//
// Created by 김재범 on 2020/11/04.
//

#ifndef ADCLASSIFIER2_KMERBUFFER_H
#define ADCLASSIFIER2_KMERBUFFER_H
#include <iostream>
#include <Kmer.h>
#include "Util.h"

class QueryKmerBuffer{
private:

public:
    QueryKmer * buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;
    QueryKmerBuffer(size_t sizeOfBuffer){
        buffer = (QueryKmer *) malloc(sizeof(QueryKmer) * sizeOfBuffer);
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    size_t reserveMemory(size_t numOfKmer){
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
        return offsetToWrite;
    };

};

class TargetKmerBuffer{
private:

public:
    TargetKmer * buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;
    TargetKmerBuffer(size_t sizeOfBuffer){
        if(sizeOfBuffer == 0){
            buffer = (TargetKmer *)malloc(sizeof(TargetKmer) * getTargetKmerBufferSize());
            bufferSize = getTargetKmerBufferSize();
            cout<<bufferSize<<endl;
        }else {
            buffer = (TargetKmer *) malloc(sizeof(TargetKmer) * sizeOfBuffer);
            bufferSize = sizeOfBuffer;
        }
        startIndexOfReserve = 0;
    };

    size_t reserveMemory(size_t numOfKmer){
        size_t offsetToWrite = __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
        return offsetToWrite;
    };

    size_t getTargetKmerBufferSize(){
        size_t memLimit = Util::getTotalSystemMemory() * 0.9;
        size_t bufferSize = memLimit / sizeof(TargetKmer);
        if(bufferSize > 10000000000){
            bufferSize = 10000000000;
        }
        return bufferSize;
    }

};
#endif //ADCLASSIFIER2_KMERBUFFER_H