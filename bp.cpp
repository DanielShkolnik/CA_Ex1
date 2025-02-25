/* 046267 Computer Architecture - Spring 2020 - HW #1 */
/* This file should hold your implementation of the predictor simulator */




#include "bp_api.h"
#include <bitset>
#include <math.h>
#include <iostream>

#define NO_SHARE 0
#define LSB_SHARE 1
#define MID_SHARE 2

//Returns the row number of the matching entry in BTB
int indx(uint32_t pc ,int btbSize){
	uint32_t indx = pc >> 2;
	indx = indx << (32 - (int)log2(btbSize));
	indx = indx >> (32 - (int)log2(btbSize));
	return indx;
}

//Returns the state_machine_number of the matching state_machine
int sharedHistory( uint32_t pc, int history, unsigned historySize, int shared, int btbSize){
	if(shared==NO_SHARE){
		return history;
	}
	else if(shared==LSB_SHARE){
		uint32_t indx = pc >> 2;
		indx = indx << (32 - historySize);
		indx = indx >> (32 - historySize);
		return indx^history;
	}
	else{
		uint32_t indx = pc >> 16;
		indx = indx << (32 - historySize);
		indx = indx >> (32 - historySize);
		return indx^history;
	}
}

//Returns the size of bits allocated in branch-predictor
unsigned calculateSize(unsigned btbSize, unsigned historySize, unsigned tagSize, bool isGlobalHist,
				  bool isGlobalTable){
	if(isGlobalHist && isGlobalTable){
		return (unsigned)(historySize+2*pow(2,historySize)+btbSize*(tagSize+30));
	}
	else if(!isGlobalHist && isGlobalTable){
		return (unsigned)(btbSize*((tagSize+30)+historySize)+2*pow(2,historySize));
	}
	else if(isGlobalHist && !isGlobalTable){
		return (unsigned)(historySize+btbSize*((tagSize+30)+2*pow(2,historySize)));
	}
	else{
		return (unsigned)(btbSize*((tagSize+30)+historySize+2*pow(2,historySize)));
	}
}

//The data held in each row in the btb
class BtbEntry{
public:
	int tag;
	uint32_t target;
	BtbEntry():tag(-1) ,target(0){}
};
//btb table
class BTB{
public:
	BtbEntry* btb;
	unsigned btbSize;
	unsigned tagSize;
	BTB(unsigned btbSize ,unsigned tagSize):btb(new BtbEntry[btbSize]) ,btbSize(btbSize) ,tagSize(tagSize){}
	~BTB(){
	    if(btb!= nullptr) delete[] btb;
	}
	void addEntry(uint32_t pc ,uint32_t targetPc){
		int i = indx(pc ,btbSize);
		btb[i].target = targetPc;

		uint32_t tagt = pc << (32 -(2 + (int)log2(this->btbSize) + this->tagSize));
        //uint32_t tagt = pc << (32 -(2 + this->tagSize));

		btb[i].tag = tagt >> (32 -(2 + (int)log2(this->btbSize) + this->tagSize) + 2 + (int)log2(this->btbSize));
        //btb[i].tag = tagt >> (32 -(2 + this->tagSize) + 2);


	}
	void updateTarget(uint32_t pc , uint32_t targetPc){
        int i = indx(pc ,btbSize);
        btb[i].target = targetPc;
	}
};

//fsm table
class FSM{
public:
	bool isGlobalTable;
	unsigned rows;
	unsigned columns;
	unsigned fsmState;
	int* fsm;
	FSM(bool isGlobalTable , unsigned btbSize , unsigned historySize , unsigned fsmState):
	isGlobalTable(isGlobalTable) ,rows(-1) ,columns(pow(2,historySize)) ,fsmState(fsmState), fsm(nullptr) {
		if(isGlobalTable){
			columns = pow(2 ,historySize);
			rows = 1;
			fsm = new int[columns];
			for(unsigned i = 0 ; i < columns ; i++)
				fsm[i] = fsmState;
		}else{
			columns = pow(2 ,historySize);
			rows = btbSize;
			fsm = new int[rows*columns];
			for(unsigned i = 0 ; i < rows ; i++){
				for(unsigned j = 0 ; j < columns ; j++){
					fsm[(i*columns)+j] = fsmState;
				}
			}
		}
	}
	~FSM(){
	    if(fsm!= nullptr) delete[] fsm;
	}
	bool isTaken(int i , unsigned history){
		if(isGlobalTable){
			if(fsm[history] > 1) {
				return true;
			}else{
				return  false;
			}
		}else{
			if(fsm[(i*columns)+history] > 1) {
				return true;
			}else{
				return  false;
			}
		}
	}
	void strengthen(int row , unsigned history){
		if(isGlobalTable){
            if(fsm[history] == 1) {
                fsm[history]--;
            }
            else if(fsm[history] == 2){
                fsm[history]++;
            }
		}else{
            if(fsm[(row*columns)+history] == 1) {
                fsm[(row*columns)+history]--;
            }
            else if(fsm[(row*columns)+history] ==2) {
                fsm[(row*columns)+history]++;
            }
		}
	}
	void weaken(int row , unsigned history){
		if(isGlobalTable){
			if(fsm[history] == 1 || fsm[history] == 0 ) {
				fsm[history]++;
			}
			else if(fsm[history] == 2 || fsm[history] == 3 ){
                fsm[history]--;
			}
		}else{
			if(fsm[(row*columns)+history] ==1 || fsm[(row*columns)+history] ==0) {
				fsm[(row*columns)+history]++;
			}
			else if(fsm[(row*columns)+history] ==2 || fsm[(row*columns)+history] ==3) {
                fsm[(row*columns)+history]--;
            }
		}
	}

	void reset(int row){
		if(!isGlobalTable) {
			for (unsigned i = 0; i < columns; i++) {
				fsm[(row * columns) + i] = fsmState;
			}
		}
	}

	void print(int row){
		if(isGlobalTable) {
			std::cout << "fsm - ";
			for(unsigned i = 0 ; i < columns ; i++)
				std::cout << fsm[i] << ",";
		}
		else{
			std::cout << "fsm - ";
			for(unsigned i = 0 ; i < columns ; i++){
				std::cout << fsm[(row*columns)+i] << ",";
			}
		}
	}

};

//BHR
class History{
public:
	unsigned* history;
	int size;
	bool isGlobalHist;
	int btbSize;
	int fsmState;
	History() = default;
	History(bool isGlobalHist , unsigned historySize ,int btbSize, int fsmState):history(nullptr) ,size(historySize) ,isGlobalHist(isGlobalHist) ,btbSize(btbSize), fsmState(fsmState){
		if(!isGlobalHist) {
			history = new unsigned[btbSize];
			for(int i = 0 ; i < btbSize ; i++)
				history[i] = 0;
		}else{
			history = new unsigned(0);
		}
	}
	~History(){
	    if(history!= nullptr){
	        if(isGlobalHist) delete history;
	        else delete[] history;
	    }
	}
	unsigned operator()(int i){
		if(isGlobalHist) {
			return *history;
		}else{
			return history[i];
		}
	}
	void reset(int row){
		if(!isGlobalHist) history[row]=0;
	}
	void print(int row){
		if(isGlobalHist) {
			std::cout << "history - " << history[0] << "  ";
		}
		else{
			std::cout << "history - " << history[row] << "  ";
		}
	}
	void update(int i ,bool taken){
		unsigned* to_update = nullptr;
		if(!isGlobalHist) {
			to_update = &history[i];
		}else{
			to_update = history;
		}
		if(size == 1) {
			std::bitset<1> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 2) {
			std::bitset<2> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}if(size == 3) {
			std::bitset<3> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 4) {
			std::bitset<4> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 5) {
			std::bitset<5> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 6) {
			std::bitset<6> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 7) {
			std::bitset<7> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 8) {
			std::bitset<8> temp(*to_update);
			temp <<= 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
	}

};

class BranchPredictor {
public:
	bool isGlobalHist;
	bool isGlobalTable;
	BTB btb;
	FSM fsm;
	History history;
	uint32_t pc;
	unsigned tagSize;
	//unsigned fsmState;
	int shared;
	unsigned btbSize;
	int historySize;

	BranchPredictor(unsigned btbSize, unsigned historySize, unsigned tagSize, unsigned fsmState,
					bool isGlobalHist, bool isGlobalTable, int shared) : isGlobalHist(isGlobalHist), isGlobalTable(isGlobalTable),
																		 btb(btbSize ,tagSize),
																		 fsm(isGlobalTable ,btbSize ,historySize ,fsmState),
																		 history(isGlobalHist ,historySize ,btbSize, fsmState),
																		 pc(0) ,tagSize(tagSize) ,shared(shared) ,btbSize(btbSize),historySize(historySize) {
	}

	bool doesExist(uint32_t pc){
		int btbRow=indx(pc,this->btbSize);


		uint32_t tag = pc << (32 -(2 + (int)log2(this->btbSize) + this->tagSize));
        //uint32_t tag = pc << (32 -(2 + this->tagSize));

		tag = tag >> (32 -(2 + (int)log2(this->btbSize) + this->tagSize) + 2 + (int)log2(this->btbSize));
        //tag = tag >> (32 -(2 + this->tagSize) + 2);

		if(this->btb.btb[btbRow].tag != -1 && (uint32_t)(this->btb.btb[btbRow].tag) == tag) return true;
		return false;

	}

	void print(uint32_t pc){
		int row = indx(pc,this->btbSize);
		this->history.print(row);
		this->fsm.print(row);
		std::cout << std::endl;
	}

};

static BranchPredictor* bp = nullptr;
static SIM_stats stats;

int BP_init(unsigned btbSize, unsigned historySize, unsigned tagSize, unsigned fsmState,
			bool isGlobalHist, bool isGlobalTable, int Shared){
	try{
		bp = new BranchPredictor(btbSize, historySize, tagSize, fsmState, isGlobalHist, isGlobalTable, Shared);
		stats.size = calculateSize(btbSize, historySize, tagSize, isGlobalHist, isGlobalTable);
		stats.flush_num = 0;
		stats.br_num = 0;
	}catch(std::bad_alloc&) {
		return -1;
	}
	return 0;
}

bool BP_predict(uint32_t pc, uint32_t *dst){
	if(bp->doesExist(pc)){
		int i = indx(pc ,bp->btbSize);
		int sharedi=sharedHistory(pc,bp->history(i),bp->historySize,bp->shared,bp->btbSize);
		bool is_taken = bp->fsm.isTaken(i ,sharedi);
		if(is_taken){
			*dst = bp->btb.btb[i].target;
		}else{
			*dst = pc + 4;
		}
		return is_taken;
	}else {
		*dst = pc + 4;
		return false;
	}
}

void BP_update(uint32_t pc, uint32_t targetPc, bool taken, uint32_t pred_dst){
	stats.br_num++;
	int i = indx(pc ,bp->btbSize);
	int sharedi=sharedHistory(pc,bp->history(i),bp->historySize,bp->shared,bp->btbSize);

	if(bp->doesExist(pc)){
		//bp->print(pc);
		bool isTaken=bp->fsm.isTaken(i ,sharedi);
		if(isTaken == taken){
			bp->fsm.strengthen(i ,sharedi);
		}else{
			bp->fsm.weaken(i ,sharedi);
		}

        if  ( (isTaken!=taken) || (taken && pred_dst!=targetPc)){
            stats.flush_num++;
        }

		bp->btb.updateTarget(pc,targetPc);

	}else {
		bp->btb.addEntry(pc, targetPc);
		//bp->print(pc);
        bp->fsm.reset(i);
		bp->history.reset(i);
        sharedi=sharedHistory(pc,bp->history(i),bp->historySize,bp->shared,bp->btbSize);
		if(bp->fsm.isTaken(i ,sharedi) == taken){
			bp->fsm.strengthen(i ,sharedi);
		}
		else {
			bp->fsm.weaken(i ,sharedi);
		}

        if (taken && pred_dst!=targetPc){
            stats.flush_num++;
        }
	}
	bp->history.update(i ,taken);
	//bp->print(pc);
}

void BP_GetStats(SIM_stats *curStats){
	*curStats = stats;
	if(bp!= nullptr) delete bp;
}


