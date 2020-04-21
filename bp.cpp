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
	uint32_t indx = pc >> 2;
	indx = indx << (32 - (int)log2(btbSize));
	indx = indx >> (32 - (int)log2(btbSize));
	return indx;
}

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


class BtbEntry{
public:
	//bool taken;
	int tag;
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

		//uint32_t tagt = pc << (32 -(2 + (int)log2(this->btbSize) + this->tagSize));
        uint32_t tagt = pc << (32 -(2 + this->tagSize));

		//btb[i].tag = tagt >> (32 -(2 + (int)log2(this->btbSize) + this->tagSize) + 2 + (int)log2(this->btbSize));
        btb[i].tag = tagt >> (32 -(2 + this->tagSize) + 2);


	}
};
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
		    //std::cout << "row: " << row << "  history: " << history << std::endl;
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
            //std::cout << "row: " << row << "  history: " << history << std::endl;
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
		//else{
		//	for(unsigned i = 0 ; i < columns ; i++){
		//		fsm[i] = fsmState;
		//	}

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
			std::cout << std::endl;
		}
	}

};
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


		//uint32_t tag = pc << (32 -(2 + (int)log2(this->btbSize) + this->tagSize));
        uint32_t tag = pc << (32 -(2 + this->tagSize));

		//tag = tag >> (32 -(2 + (int)log2(this->btbSize) + this->tagSize) + 2 + (int)log2(this->btbSize));
        tag = tag >> (32 -(2 + this->tagSize) + 2);

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
		//std::cout << "  does exist  " << std::endl;
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
    //std::cout <<"row: " << i << std::endl;
    //std::cout <<"bp->history(i): " << bp->history(i) << std::endl;
	int sharedi=sharedHistory(pc,bp->history(i),bp->historySize,bp->shared,bp->btbSize);
    //std::cout <<"bp->history(i): " << sharedi << std::endl;
	if(bp->doesExist(pc) && pred_dst==targetPc){
		//bp->print(pc);
		if(bp->fsm.isTaken(i ,sharedi) == taken){
			bp->fsm.strengthen(i ,sharedi);
		}else{
			//std::cout << pc << "  flushed  " << std::endl;
			stats.flush_num++;
			bp->fsm.weaken(i ,sharedi);
		}
	}else {
		//if(!bp->isGlobalTable || bp->btb.btb[i].tag!=-1) bp->fsm.reset(i);
		bp->btb.addEntry(pc, targetPc);
        bp->fsm.reset(i);
		bp->history.reset(i);
		//bp->print(pc);
        sharedi=sharedHistory(pc,bp->history(i),bp->historySize,bp->shared,bp->btbSize);
		if(bp->fsm.isTaken(i ,sharedi) == taken){
			bp->fsm.strengthen(i ,sharedi);
		}
		else {
			//std::cout << pc << "  flushed  " << std::endl;
			stats.flush_num++;
			bp->fsm.weaken(i ,sharedi);
		}
	}
	bp->history.update(i ,taken);
	//bp->print(pc);
}

void BP_GetStats(SIM_stats *curStats){

	*curStats = stats;
}





/*

#include "bp_api.h"


#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>

typedef enum {LOCAL=0,GLOBAL=1} table_type;
typedef  enum {NOT_USING_SHARE=0,USING_SHARE_LSB=1,USING_SHARE_MID=2} share_type;
typedef  enum {SNT,WNT,WT,ST} FSM_STATE;



uint32_t cut_address(uint32_t adress , uint32_t chunk ,uint32_t div){
    assert(chunk<32);
    uint32_t res =adress;
    res = res >> div;
    uint32_t mask = 0xFFFFFFFF;
    assert(32-chunk>=0);
    mask = mask >> (32-chunk);
    return mask & res ;
}

unsigned hash_func(uint32_t cur_pc, unsigned history , unsigned history_size ,share_type share_type ){
    unsigned res;
    unsigned mask=1;
    for (unsigned i = 0; i < history_size-1; ++i) {
        mask*=2;
        mask++;
    }
    unsigned div;
    switch (share_type){
        case  NOT_USING_SHARE :
            return history;
        case USING_SHARE_LSB:
            div=2;
            break;
        case USING_SHARE_MID:
            div=16;
            break;
    }
    res=cut_address(cur_pc,history_size,div);
    return (res & mask) ^ history;
}



// fsm class
class fsm
{
public:
    FSM_STATE cur;
    FSM_STATE default_state;

public:
    fsm (FSM_STATE default_s): cur(default_s) ,default_state(default_s) {}
    void update(bool input){
        switch (cur) {
            case ST :
                cur = (input) ? ST : WT;
                break;
            case WT :
                cur = (input) ? ST : WNT;
                break;
            case WNT :
                cur = (input) ? WT : SNT;
                break;
            case SNT :
                cur = (input) ? WNT : SNT;
                break;
        }
    }
    FSM_STATE get(){
        return cur;
    }
    void reset (){
        cur=default_state;
    }

};


class row
{
    uint32_t tag = 1;
    std::vector<fsm>* fsm_table_pointer=nullptr;
    uint32_t target_pc;
    unsigned* history =nullptr;
    table_type history_rule;
    table_type fsm_rule;
    unsigned history_size;
    FSM_STATE default_state;

public:
    bool valid_bit=false;

    row(uint32_t target_pc,unsigned  history_size,table_type history_rule,table_type fsm_rule,std::vector<fsm>* fsm_table_pointer
            ,unsigned* history,FSM_STATE default_s): fsm_table_pointer(fsm_table_pointer) ,target_pc(target_pc)
            ,history(history),history_rule(history_rule) , fsm_rule(fsm_rule),history_size(history_size),default_state(default_s)

    {

        if(history_rule == LOCAL){
            //history = new unsigned(0);        //local history
        } else{
            this->history = history; //global history
        }
        if(fsm_rule == LOCAL){
            //fsm_table_pointer = new std::vector<fsm>(std::pow(2,history_size),fsm(default_s)); //local fsm
        } else{
            this->fsm_table_pointer =fsm_table_pointer; //global fsm
        }
    }
    ~row(){
        if(fsm_rule == LOCAL){
            if(fsm_table_pointer !=nullptr	){
                delete  fsm_table_pointer;
            }
        }
        if (history_rule == LOCAL){
            if(history != nullptr){
                delete  history;
            }
        }
    }
    row(const row& r)
    {
        tag=r.tag;
        default_state=r.default_state;
        target_pc=r.target_pc;
        history_rule=r.history_rule;
        fsm_rule=r.fsm_rule;
        history_size=r.history_size;
        if(history_rule == LOCAL){
            history = new unsigned(0);        //local history
        } else{
            history = r.history; //global history
        }
        if (fsm_rule == LOCAL){
            fsm_table_pointer = new std::vector<fsm>(std::pow(2,history_size),fsm(r.default_state));
        }
        else{
            fsm_table_pointer =r.fsm_table_pointer;
        }
    }
    row& operator=(const row& r){
        if (this == &r){
            return *this;
        }
        delete history;
        delete fsm_table_pointer;
        return *this;
    }
    uint32_t get_tag(){
        return tag;
    }

    uint32_t get_dst(){
        return target_pc;
    }

    void reset_history (){
        *history=0;
    }

    void update_history (bool input){
        *history=*(history)*2+input;
        unsigned mask=1;
        for (unsigned i = 0; i < history_size-1; ++i) {
            mask*=2;
            mask++;
        }
        *history=*history & mask;
    }
    //assume valid tag in the input
    void update_tag(uint32_t new_tag){
        tag = new_tag;
    }
    void update_dst(uint32_t new_dst){
        target_pc = new_dst;
    }

    bool is_taken (uint32_t cur_pc,share_type share_type){
        FSM_STATE res=(fsm_table_pointer)->operator[](hash_func(cur_pc,*history,history_size,share_type)).get();
        return (res <= WNT) ? false : true ;
    }

    void reset_fsm_vector(){
        for (unsigned i=0 ; i < fsm_table_pointer->size(); i++)
        {
            fsm_table_pointer->operator[](i).reset();
        }
    }

    std::vector<fsm>* get_fsm_vector(){
        return fsm_table_pointer;
    }

    unsigned get_history(){
        return *history;
    }
};




class btb
{
public:

    std::vector<row>* b_data =nullptr;
    unsigned  history_size;
    unsigned  tag_size;
    FSM_STATE default_state;
    table_type history_rule;    //global or local gistory table
    table_type fsm_rule;        //global or local fsm
    share_type share_t;
    std::vector<fsm>* fsm_g = nullptr;
    unsigned* history_g = nullptr;
    SIM_stats stats;


public:
    btb (){}
    btb (unsigned btbSize, unsigned historySize, unsigned tagSize, unsigned fsmState,
         bool isGlobalHist, bool isGlobalTable, int Shared):
            history_size(historySize),tag_size(tagSize),default_state((FSM_STATE)fsmState),history_rule((table_type)isGlobalHist),
            fsm_rule((table_type)isGlobalTable),share_t((share_type)Shared)	 {

        if(fsm_rule == GLOBAL){
            fsm_g = new std::vector<fsm>(std::pow(2,historySize),fsm((FSM_STATE)fsmState));
        }
        if(history_rule == GLOBAL){
            history_g = new unsigned(0); //global history
        }

        b_data=new std::vector<row>(btbSize,row(0,historySize,history_rule,fsm_rule,fsm_g,history_g,default_state));
        stats.br_num=0;
        stats.flush_num=0;
        stats.size=btb_size_calc();

    }
    ~btb(){
        if(fsm_rule == GLOBAL){
            delete  fsm_g;
        }

        if(history_rule == GLOBAL){
            delete history_g;
        }
        delete b_data;
    }
    SIM_stats get_stats(){
        return stats;
    }

    //get the prediction address and the Taken/Not Taken prediction (Decode stage)
    bool get_dst(uint32_t pc , uint32_t* dst){
        row& req_row = b_data->operator [](cut_address(pc,log2((double)b_data->size()),2)); //div the offset 2 times.
        if(!req_row.valid_bit || ( req_row.get_tag() != cut_address(pc,tag_size,2))) {
            *dst = pc+4;
            return false;
        }
        if (! req_row.is_taken(pc,share_t)){
            *dst = pc+4;
            return false;
        }
        *dst=req_row.get_dst();
        return true;
    }
    //update the btb (EXE stage)
    void update(uint32_t pc, uint32_t targetPc, bool taken, uint32_t pred_dst) {
        stats.br_num++;
        if  ( (taken && pred_dst==pc+4) || (!taken && pred_dst!=pc+4) || (taken && pred_dst!=targetPc)){
            stats.flush_num++;
        }
        row& row = b_data->operator [](cut_address(pc,log2((double)b_data->size()),2)); //div the offset 2 times.
        if(!known_branch(pc,targetPc,row)){

            row.update_tag(cut_address(pc,tag_size,2));
            if(history_rule == LOCAL){
                row.reset_history();
            }
            if(fsm_rule == LOCAL){
                row.reset_fsm_vector();
            }
        }
        //step -update the current line to the btb
        row.update_dst(targetPc);
        row.get_fsm_vector()->operator[](hash_func(pc,row.get_history(),history_size,share_t)).update(taken); //update FSM
        row.update_history(taken); //update history
        row.valid_bit=true;
    }



    bool known_branch(uint32_t pc, uint32_t targetPc ,row row){
        if(row.get_tag() != cut_address(pc,tag_size,2)){
            return false; // tag conflict
        }
        return true;
    }

    //calc the size of the btb
    int btb_size_calc(){
        int btb_size=b_data->size();
        int tag_table_size=(btb_size)*(tag_size + 30);
        int history_table_size=0;
        int fsm_table_size=0;
        switch (history_rule){
            case LOCAL :
                history_table_size=btb_size*history_size;
                break;
            case GLOBAL :
                history_table_size=history_size;
                break;
        }
        switch (fsm_rule){
            case LOCAL :
                fsm_table_size=btb_size* (std::pow(2,history_size+1));
                break;
            case GLOBAL :
                fsm_table_size= (std::pow(2,history_size+1));
                break;
        }
        return history_table_size+fsm_table_size+tag_table_size;
    }
};

btb* global_var_btb= nullptr; //global btb


int BP_init(unsigned btbSize, unsigned historySize, unsigned tagSize, unsigned fsmState,
            bool isGlobalHist, bool isGlobalTable, int Shared){
    global_var_btb = new  btb(btbSize,historySize,tagSize,fsmState,isGlobalHist,isGlobalTable,Shared);
    if(global_var_btb == nullptr){
        return -1; //memory error
    }
    return 0;
}


bool BP_predict(uint32_t pc, uint32_t *dst){
    return global_var_btb->get_dst(pc,dst);
}

void BP_update(uint32_t pc, uint32_t targetPc, bool taken, uint32_t pred_dst){
    std::cout << "history - " << global_var_btb->history_g[0] << "  ";
    std::cout << "fsm - ";
    for(auto i=global_var_btb->fsm_g->begin(); i!=global_var_btb->fsm_g->end(); i++){
        std::cout << i->cur << ",";
    }
    std::cout << std::endl;
    global_var_btb->update(pc,targetPc,taken,pred_dst);
    std::cout << "history - " << global_var_btb->history_g[0] << "  ";
    std::cout << "fsm - ";
    for(auto i=global_var_btb->fsm_g->begin(); i!=global_var_btb->fsm_g->end(); i++){
        std::cout << i->cur << ",";
    }
    std::cout << std::endl;
    return ;
}

void BP_GetStats(SIM_stats *curStats){
    *curStats = global_var_btb->get_stats();
    delete  global_var_btb;
    global_var_btb= nullptr;
    return;
}
*/



