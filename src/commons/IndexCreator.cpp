//
// Created by KJB on 01/09/2020.
//

#include "IndexCreator.h"
//#include "NcbiTaxonomy.h"
IndexCreator::IndexCreator()
{
    kmerExtractor = new KmerExtractor();
}

IndexCreator::~IndexCreator() { delete kmerExtractor;}
void IndexCreator::startIndexCreating2(ifstream & targetFile, const char * outputFileName, vector<int> & taxIdList)
{
    string buffer;
    string forwardRead;
    string reverseComplimentRead;
    string reads[2];
    ExtractStartPoint ESP = {0 ,0};
    size_t bufferIdx = 0;

    struct MmapedData<char> seqFile = mmapData<char>( "/Users/kjb/Desktop/ADclassifier/tengenome/tengenome.fna");
    size_t maxNuc = seqFile.fileSize/sizeof(char);
    vector<SeqSegment> seqSegments;
    size_t start = 0;
    size_t end = 0;
    SeqSegment temp;
    for(size_t i = 0; i < maxNuc; i++)
    {
        if(seqFile.data[i] == '>')
        {
            end = i - 2;
            temp = {start, end};
            seqSegments.push_back(temp); // the first push_back is a garbage.
            while(seqFile.data[i] != '\n')
            {
                cout<<seqFile.data[i];
                i++;
            }
            cout<<endl;
            start = i + 1;
        }
    }
    temp = {start, maxNuc - 2};
    seqSegments.push_back(temp);
    cout<<seqSegments.size()<<endl;
    cout<<seqSegments[18].start<<" "<<seqSegments[18].end<<endl;
    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);

    int seqID = 0;
    for(size_t i = 1 ; i < seqSegments.size(); i++)
    {
        kmerExtractor -> dna2aa2(seqSegments[i], seqFile);
        ESP = kmerExtractor->fillKmerBuffer2(seqSegments[i], seqFile, kmerBuffer, seqID, bufferIdx, ESP);
        while (ESP.startOfFrame + ESP.frame != 0)
        {
            writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdList);
            ESP = kmerExtractor->fillKmerBuffer2(seqSegments[i], seqFile, kmerBuffer, seqID, bufferIdx, ESP);
        }
        seqID ++;
    }

    //flush last buffer
    writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdList);
    targetFile.close();
    free(kmerBuffer);
    munmap(seqFile.data, seqFile.fileSize + 1);
}
void IndexCreator::startIndexCreating(ifstream & targetFile, char * outputFileName, vector<int> & taxIdList)
{
    string buffer;
    string forwardRead;
    string reverseComplimentRead;
    string reads[2];
    ExtractStartPoint ESP = {0 ,0};
    size_t bufferIdx = 0;

    Kmer * kmerBuffer = (Kmer *)malloc(sizeof(Kmer) * kmerBufSize);
    int seqID = 0;
    getline(targetFile, buffer);
    while(targetFile)
    {
        getline(targetFile, buffer);
        if(buffer[0] == '>'){
            reverseComplimentRead = kmerExtractor->reverseCompliment(forwardRead);
            reads[0] = forwardRead; reads[1] = reverseComplimentRead;
            kmerExtractor->dna2aa(forwardRead, reverseComplimentRead);
//            seqID = taxIdList[idIdx]; idIdx++;
            ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
            while (ESP.startOfFrame + ESP.frame != 0)
            {
                writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdList);
                ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
            }
            forwardRead.clear();
            seqID ++;

            continue;
        }
        forwardRead.append(buffer);
    }
    // For last one
    reverseComplimentRead = kmerExtractor->reverseCompliment(forwardRead);
    reads[0] = forwardRead; reads[1] = reverseComplimentRead;
    kmerExtractor->dna2aa(forwardRead, reverseComplimentRead);
    //seqID = taxIdList[idIdx]; idIdx++;
    ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
    while (ESP.startOfFrame + ESP.frame != 0)
    {
        writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdList);
        ESP = kmerExtractor->fillKmerBuffer(reads, kmerBuffer, seqID, bufferIdx, ESP);
    }

    //flush last buffer
    writeTargetFiles(kmerBuffer, bufferIdx, outputFileName, taxIdList);
    targetFile.close();
    free(kmerBuffer);

}

void IndexCreator::writeTargetFiles(Kmer *kmerBuffer, size_t & bufferIdx, const char * outputFileName, vector<int> & taxIdList)
{
    char suffixedDiffIdxFileName[100];
    char suffixedInfoFileName[100];
    sprintf(suffixedDiffIdxFileName,"%s_diffIdx_%zu", outputFileName,numOfFlush);
    sprintf(suffixedInfoFileName,"%s_info_%zu", outputFileName,numOfFlush);
    numOfFlush++;
    FILE * diffIdxFile = fopen(suffixedDiffIdxFileName, "wb");
    FILE * idAndPosFile = fopen(suffixedInfoFileName, "wb");
    cout<<suffixedInfoFileName<<endl;

    uint16_t *kmerLocalBuf = (uint16_t *)malloc(sizeof(uint16_t) * kmerBufSize);
    size_t localBufIdx = 0;

    if (diffIdxFile == NULL || idAndPosFile == NULL){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    uint64_t lastKmer = 0;

    sort(kmerBuffer, kmerBuffer + bufferIdx, [](const Kmer & a, const Kmer & b) {
        return a.ADkmer < b.ADkmer || (a.ADkmer == b.ADkmer && a.info.sequenceID < b.info.sequenceID);
    });
    Kmer lookingKmer = kmerBuffer[0];
    size_t write = 0;
    int endFlag = 0;


    for(size_t i = 1 ; i < bufferIdx ; i++) {
        while(taxIdList[lookingKmer.info.sequenceID] == taxIdList[kmerBuffer[i].info.sequenceID]){
            if (lookingKmer.ADkmer != kmerBuffer[i].ADkmer) {
                break;
            }
            i++;
            if(i == bufferIdx)
            {
                endFlag = 1;
                break;
            }
        }

        fwrite(&lookingKmer.info, sizeof(KmerInfo), 1, idAndPosFile);
        write++;
        writeKmerDiff(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);

        if(endFlag == 1) break;

        lastKmer = lookingKmer.ADkmer;
        lookingKmer = kmerBuffer[i];
    }
    if(!((kmerBuffer[bufferIdx - 2].ADkmer == kmerBuffer[bufferIdx - 1].ADkmer) &&
        (kmerBuffer[bufferIdx - 2].info.sequenceID == kmerBuffer[bufferIdx - 1].info.sequenceID))){
        fwrite(&lookingKmer.info, sizeof(KmerInfo), 1, idAndPosFile);
        write++;
        writeKmerDiff(lastKmer, lookingKmer.ADkmer, diffIdxFile, kmerLocalBuf, localBufIdx);
    }

    cout<<"total k-mer count: "<<bufferIdx<<endl;
    cout<<"written k-mer count: "<<write<<endl;
    flushKmerBuf(kmerLocalBuf, diffIdxFile, localBufIdx);

    free(kmerLocalBuf);
    fclose(diffIdxFile);
    fclose(idAndPosFile);
    bufferIdx = 0;
}

void IndexCreator::writeKmerDiff(uint64_t lastKmer, uint64_t & entryToWrite, FILE* handleKmerTable, uint16_t *kmerBuf, size_t & localBufIdx ){
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[5];
    int idx = 3;

    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;

    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }

    writeKmer(kmerBuf, handleKmerTable, (buffer + idx + 1), (4 - idx), localBufIdx);
}

void IndexCreator::flushKmerBuf(uint16_t *buffer, FILE *handleKmerTable, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(uint16_t), localBufIdx, handleKmerTable);
    localBufIdx = 0;
}

void IndexCreator::writeKmer(uint16_t *buffer, FILE* handleKmerTable, uint16_t *toWrite, size_t size, size_t & localBufIdx ) {
    if (localBufIdx + size >= kmerBufSize) {
        flushKmerBuf(buffer, handleKmerTable, localBufIdx);
    }
    memcpy(buffer + localBufIdx, toWrite, sizeof(uint16_t) * size);
    localBufIdx += size;
}

void IndexCreator::writeInfo(KmerInfo * entryToWrite, FILE * infoFile, KmerInfo * infoBuffer, size_t & infoBufferIdx)
{
    if (infoBufferIdx >= kmerBufSize) {
        flushInfoBuf(infoBuffer, infoFile, infoBufferIdx);
    }
    memcpy(infoBuffer + infoBufferIdx, entryToWrite, sizeof(KmerInfo));
    infoBufferIdx++;
}
void IndexCreator::flushInfoBuf(KmerInfo * buffer, FILE * infoFile, size_t & localBufIdx ) {
    fwrite(buffer, sizeof(KmerInfo), localBufIdx, infoFile);
    localBufIdx = 0;
}
int IndexCreator::getNumOfFlush()
{
    return numOfFlush;
}