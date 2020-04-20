/* 046267 Computer Architecture - Spring 2020 - HW #1 */
/* This file should hold your implementation of the predictor simulator */

#include "bp_api.h"
#include <bitset>
#include <math.h>
#include <iostream>

#define NO_SHARE 0
#define LSB_SHARE 1
#define MID_SHARE 2

int indx(uint32_t pc ,int btbSize){
	uint32_t indx = pc << 2;
	indx = indx << (32 - (int)log2(btbSize));
	indx = indx >> (32 - (int)log2(btbSize));
	return indx;
}

unsigned calculateSize(unsigned btbSize, unsigned historySize, unsigned tagSize, bool isGlobalHist,
				  bool isGlobalTable){
	if(isGlobalHist && isGlobalTable){
		return (unsigned)(historySize+2*pow(2,historySize)+btbSize*tagSize);
	}
	else if(!isGlobalHist && isGlobalTable){
		return (unsigned)(btbSize*(tagSize+historySize)+2*pow(2,historySize));
	}
	else if(isGlobalHist && !isGlobalTable){
		return (unsigned)(historySize+btbSize*(tagSize+2*pow(2,historySize)));
	}
	else{
		return (unsigned)(btbSize*(tagSize+historySize+2*pow(2,historySize)));
	}
}


class BtbEntry{
public:
	//bool taken;
	uint32_t tag;
	uint32_t target;
	//int history;
	//int* states;
	BtbEntry():tag(-1) ,target(0){}
};
class BTB{
public:
	BtbEntry* btb;
	unsigned btbSize;
	unsigned tagSize;
	BTB(unsigned btbSize ,unsigned tagSize):btb(new BtbEntry[btbSize]) ,btbSize(btbSize) ,tagSize(tagSize){}
	void addEntry(uint32_t pc ,uint32_t targetPc){
		int i = indx(pc ,btbSize);
		btb[i].target = targetPc;
		uint32_t tagt = pc << (32 -(2 + (int)log2(this->btbSize) + this->tagSize));
		std::cout << tagt << std::endl;
		btb[i].tag = tagt >> (32 -(2 + (int)log2(this->btbSize) + this->tagSize) + 2 + (int)log2(this->btbSize));
		std::cout << btb[i].tag << std::endl;

	}
};
class FSM{
public:
	bool isGlobalTable;
	unsigned column_size;
	unsigned row_size;
	unsigned fsmState;
	int* fsm;
	FSM(bool isGlobalTable , unsigned btbSize , unsigned historySize , unsigned fsmState):isGlobalTable(isGlobalTable) ,column_size(-1) ,row_size(pow(2,historySize)) ,fsm(nullptr) ,fsmState(fsmState){
		if(isGlobalTable){
			row_size = pow(2 ,historySize);
			column_size = 1;
			fsm = new int[row_size];
			for(int i = 0 ; i < row_size ; i++)
				fsm[i] = fsmState;
		}else{
			row_size = pow(2 ,historySize);
			column_size = btbSize;
			fsm = new int[column_size*row_size];
			for(int i = 0 ; i < column_size ; i++){
				for(int j = i ; j < row_size ; j++){
					fsm[(i*column_size)+j] = fsmState;
				}
			}
		}
	}
	bool isTaken(int i , unsigned history){
		if(isGlobalTable){
			if(fsm[history] > 1) {
				return true;
			}else{
				return  false;
			}
		}else{
			if(fsm[(i*column_size)+history] > 1) {
				return true;
			}else{
				return  false;
			}
		}
	}
	void strengthen(int i , unsigned history){
		if(isGlobalTable){
			if(fsm[history] < 3) {
				fsm[history]++;
			}
		}else{
			if(fsm[(i*column_size)+history] < 3) {
				fsm[(i*column_size)+history]++;
			}
		}
	}
	void weaken(int i , unsigned history){
		if(isGlobalTable){
			if(fsm[history] > 0) {
				fsm[history]--;
			}
		}else{
			if(fsm[(i*column_size)+history] > 0) {
				fsm[(i*column_size)+history]--;
			}
		}
	}

};
class History{
public:
	unsigned* history;
	int size;
	bool isGlobalHist;
	int btbSize;
	History() = default;
	History(bool isGlobalHist , unsigned historySize ,int btbSize):history(nullptr) ,size(historySize) ,isGlobalHist(isGlobalHist) ,btbSize(btbSize){
		if(!isGlobalHist) {
			history = new unsigned[btbSize];
			for(int i = 0 ; i < btbSize ; i++)
				history[i] = 0;
		}else{
			history = new unsigned(0);
		}
	}
	unsigned operator()(int i){
		if(isGlobalHist) {
			return *history;
		}else{
			return history[i];
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
			temp << 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 2) {
			std::bitset<2> temp(*to_update);
			temp << 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}if(size == 3) {
			std::bitset<3> temp(*to_update);
			temp << 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 4) {
			std::bitset<4> temp(*to_update);
			temp << 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 5) {
			std::bitset<5> temp(*to_update);
			temp << 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 6) {
			std::bitset<6> temp(*to_update);
			temp << 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 7) {
			std::bitset<7> temp(*to_update);
			temp << 1;
			temp[0] = taken;
			*to_update = temp.to_ulong();
			return;
		}
		if(size == 8) {
			std::bitset<8> temp(*to_update);
			temp << 1;
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

	BranchPredictor(unsigned btbSize, unsigned historySize, unsigned tagSize, unsigned fsmState,
					bool isGlobalHist, bool isGlobalTable, int Shared) : isGlobalHist(isGlobalHist), isGlobalTable(isGlobalTable),
																		 btb(btbSize ,tagSize),
																		 fsm(isGlobalTable ,btbSize ,historySize ,fsmState),
																		 history(isGlobalHist ,historySize ,btbSize),
																		 pc(0) ,tagSize(tagSize) ,shared(shared) ,btbSize(btbSize) {
	}

	bool doesExist(uint32_t pc){
		int btbRow=indx(pc,this->btbSize);
		uint32_t tag = pc << (32 -(2 + (int)log2(this->btbSize) + this->tagSize));
		tag = tag >> (32 -(2 + (int)log2(this->btbSize) + this->tagSize) + 2 + (int)log2(this->btbSize));
		if(this->btb.btb[btbRow].tag == tag) return true;
		return false;

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
		std::cout << "  does exist  " << std::endl;
		int i = indx(pc ,bp->btbSize);
		bool is_taken = bp->fsm.isTaken(i ,bp->history(i));
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
	if(bp->doesExist(pc)){
		int i = indx(pc ,bp->btbSize);
		if(bp->fsm.isTaken(i ,bp->history(i)) == taken && targetPc == pred_dst){
			bp->fsm.strengthen(i ,bp->history(i));
		}else{
			stats.flush_num++;
			bp->fsm.weaken(i ,bp->history(i));
		}
	}else {
		bp->btb.addEntry(pc, targetPc);
	}
	int i = indx(pc ,bp->btbSize);
	bp->history.update(i ,taken);

}

void BP_GetStats(SIM_stats *curStats){

	*curStats = stats;
}

