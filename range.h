//
//  range.h
//  Savitzky_Golay_3D
//
//  Created by yoshi on 03/02/2023.
//

#ifndef range_h
#define range_h

class Range
{
public:
    Range(int xsize);
    Range(){}
    int init(int size){
        size_ = size;
        return 0;
    }
    size_t operator()(int i)
    {
        while(i < 0) i+= size_;
        while(i > (size_-1)) i -= size_;
        return i ;
    }
protected:
    int size_;
};

#endif /* range_h */
