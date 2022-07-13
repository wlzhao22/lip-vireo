#ifndef NNITEM_H
#define NNITEM_H


class NNItem
{
public:
	NNItem(const unsigned int index, const float dist)
	{
	    this->index = index;
	    this->val = dist;
	}
	unsigned int index;
	float val;
	float base;
	float range;
	~NNItem()
	{
	}

	//sort in ascending order
    static int comparer(const NNItem *a,const NNItem *b)
    {
         return (a->val < b->val);
    }

    //sort in descending order
    static int lgcomparer(const NNItem *a,const NNItem *b)
    {
         return (a->val > b->val);
    }
};

/********for building a heap************/
class HItem
{
public:
	HItem()
	{
	    this->val     = 0;
	    this->visited = false;
	}
	bool visited;
	unsigned int bidx;
	float val;
	~HItem()
	{
	}

	//sort in ascending order
    static int llcomp(const HItem *a,const HItem *b)
    {
         return (a->val < b->val);
    }

    //sort in descending order
    static int lgcomp(const HItem *a,const HItem *b)
    {
         return (a->val > b->val);
    }
};

#endif
