#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace sdsl;
using namespace std;

typedef sdsl::csa_wt<sdsl::wt_pc<sdsl::balanced_shape, sdsl::int_vector<(unsigned char)1>, sdsl::rank_support_v<(unsigned char)1, (unsigned char)1>, sdsl::select_support_mcl<(unsigned char)1, (unsigned char)1>, sdsl::select_support_mcl<(unsigned char)0, (unsigned char)1>, sdsl::byte_tree<false> >, 32u, 32u, sdsl::sa_order_sa_sampling<(unsigned char)0>, sdsl::isa_sampling<(unsigned char)0>, sdsl::byte_alphabet> t_csa;

size_t countWithOneError(t_csa &csa, t_csa &csaReverse, string::iterator begin,
                         string::iterator end,
                         std::vector<char> &alphabet)
{
   
    size_t x = ceil(double(end-begin)/2);
    size_t count = 0, size = 0;
    size_t leftRes, rightRes, leftNew , rightNew,leftSave=0,rightSave=0;
    string::iterator endCaseA;
    size = backward_search(csa, 0, csa.size() - 1, begin + x, end, leftRes, rightRes); //**
    endCaseA= begin+x;
    if (size != 0)
    {
        while (begin < endCaseA)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = backward_search(csa, leftRes, rightRes, alphabet[i],leftNew , rightNew);
                if (size != 0)
                {
                    if (*(endCaseA - 1) != alphabet[i])
                        count = count + backward_search(csa, leftNew, rightNew, begin, endCaseA - 1, leftNew, rightNew);
                    else{
                        leftSave=leftNew;
                        rightSave=rightNew;
                    }   
                }
            }
            endCaseA--;
            leftRes=leftSave;
            rightRes=rightSave;
        }
    }
    size_t leftB , rightB , leftRevB , rightRevB, leftBRes , rightBRes , leftRevBRes , rightRevBRes,leftBSave=0 , rightBSave=0 , leftRevBSave=0 , rightRevBSave=0;
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin, begin + x, leftB, rightB, leftRevB, rightRevB);
    string::iterator beginCaseB = begin+x ; 
    if (size != 0)
    {
        while (beginCaseB < end)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csaReverse, leftRevB, rightRevB, leftB, rightB, alphabet[i], leftRevBRes, rightRevBRes, leftBRes, rightBRes);
                if (size != 0)
                {
                    if (*(beginCaseB) != alphabet[i])
                        count = count + bidirectional_search_forward(csa, csaReverse, leftBRes, rightBRes, leftRevBRes, rightRevBRes, beginCaseB+1, end, leftBRes, rightBRes, leftRevBRes, rightRevBRes);
                    else{
                        leftRevBSave=leftRevBRes;
                        rightRevBSave=rightRevBRes;
                        leftBSave=leftBRes;
                        rightBSave=rightBRes;
                    }    
                }
            }
            beginCaseB++;
            leftRevB=leftRevBSave;
            rightRevB=rightRevBSave;
            leftB=leftBSave;
            rightB=rightBSave;

        }
    }
    return count;
}

int_vector<64> locateWithOneError(const t_csa &csa, const t_csa &csaReverse, string::iterator begin,
                                                  string::iterator end,
                                                  size_t max_locate,
                                                  std::vector<char> &alphabet)
{
    size_t x = ceil(double(end - begin) / 2);
    size_t size = 0, nextIndex = 0;
    int_vector<64> results;
    results.resize(max_locate);
    size_t leftRes, rightRes, leftNew, rightNew, leftSave=0, rightSave=0;
    string::iterator endCaseA;
    size = backward_search(csa, 0, csa.size() - 1, begin + x, end, leftRes, rightRes); //**
    endCaseA = begin + x;
    if (size != 0)
    {
        while (begin < endCaseA)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = backward_search(csa, leftRes, rightRes, alphabet[i], leftNew, rightNew);
                if (size != 0)
                {
                    if (*(endCaseA - 1) != alphabet[i])
                    {
                        size = backward_search(csa, leftNew, rightNew, begin, endCaseA - 1, leftNew, rightNew);
                        for (size_t k = 0; k < size; k++)
                            results[nextIndex + k] = csa[leftNew + k];

                        nextIndex = nextIndex + size;
                    }
                    else
                    {
                        leftSave = leftNew;
                        rightSave = rightNew;
                    }
                }
            }
            endCaseA--;
            leftRes = leftSave;
            rightRes = rightSave;
        }
    }
    size_t leftB, rightB, leftRevB, rightRevB, leftBRes, rightBRes, leftRevBRes, rightRevBRes, leftBSave=0, rightBSave=0, leftRevBSave=0, rightRevBSave=0;
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin, begin + x, leftB, rightB, leftRevB, rightRevB);
    string::iterator beginCaseB = begin + x;
    if (size != 0)
    {
        while (beginCaseB < end)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csaReverse, leftRevB, rightRevB, leftB, rightB, alphabet[i], leftRevBRes, rightRevBRes, leftBRes, rightBRes);
                if (size != 0)
                {
                    if (*(beginCaseB) != alphabet[i])
                    {
                        size = bidirectional_search_forward(csa, csaReverse, leftBRes, rightBRes, leftRevBRes, rightRevBRes, beginCaseB + 1, end, leftBRes, rightBRes, leftRevBRes, rightRevBRes);
                        for (size_t k = 0; k < size; k++)
                            results[nextIndex + k] = csa[leftBRes + k];

                        nextIndex = nextIndex + size;
                    }
                    else
                    {
                        leftRevBSave = leftRevBRes;
                        rightRevBSave = rightRevBRes;
                        leftBSave = leftBRes;
                        rightBSave = rightBRes;
                    }
                }
            }
            beginCaseB++;
            leftRevB = leftRevBSave;
            rightRevB = rightRevBSave;
            leftB = leftBSave;
            rightB = rightBSave;
        }
    }
    return results;
}




size_t TwomissCaseASubCase(
    t_csa& csa,
    size_t left,
    size_t right,
    string::iterator begin,
    string::iterator end,
    const std::vector<char>& alphabet
){
    size_t leftNew,rightNew,leftSave=0,rightSave=0;
    size_t size=0;
    size_t count=0;
    size_t leftRes=left;
    size_t rightRes=right;
    while (begin < end)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = backward_search(csa, leftRes, rightRes, alphabet[i],leftNew , rightNew);
                if (size != 0)
                {
                    if (*(end - 1) != alphabet[i])
                        count = count + backward_search(csa, leftNew, rightNew, begin, end - 1, leftNew, rightNew);
                    else{
                        leftSave=leftNew;
                        rightSave=rightNew;
                    }   
                }
            }
            end--;
            leftRes=leftSave;
            rightRes=rightSave;
        }

        return count;
}

size_t TwomissCaseBSubCase(
    t_csa& csa,
    t_csa& csaReverse,
    size_t leftB,
    size_t rightB,
    size_t leftRevB,
    size_t rightRevB,
    string::iterator begin,
    string::iterator end,
    const std::vector<char>& alphabet
)
{
    size_t leftRevBRes,rightRevBRes,leftBRes,rightBRes,leftRevBSave=0,rightRevBSave=0,leftBSave=0,rightBSave=0;
    size_t size=0;
    size_t count=0;
    size_t left=leftB;
    size_t right=rightB;
    size_t leftRev=leftRevB;
    size_t rightRev=rightRevB;
     while (begin < end)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csaReverse, leftRev, rightRev, left, right, alphabet[i], leftRevBRes, rightRevBRes, leftBRes, rightBRes);
                if (size != 0)
                {
                    if (*(begin) != alphabet[i])
                        count = count + bidirectional_search_forward(csa, csaReverse, leftBRes, rightBRes, leftRevBRes, rightRevBRes, begin+1, end, leftBRes, rightBRes, leftRevBRes, rightRevBRes);
                    else{
                        leftRevBSave=leftRevBRes;
                        rightRevBSave=rightRevBRes;
                        leftBSave=leftBRes;
                        rightBSave=rightBRes;
                    }    
                }
            }
            begin++;
            leftRev=leftRevBSave;
            rightRev=rightRevBSave;
            left=leftBSave;
            right=rightBSave;

        }

        return count;
}

size_t countWithTwoError(t_csa &csa, t_csa &csaReverse, string::iterator begin,
                         string::iterator end,
                         std::vector<char> &alphabet)
{
    size_t left;
    size_t right;
    size_t leftRev;
    size_t rightRev;
    size_t s1 = floor((double)(end-begin)/3);
    size_t s2 = end - begin - s1;
    size_t count = 0, size = 0;
    size_t leftRes, rightRes, leftRevRes, rightRevRes, leftSave=0,rightSave=0,leftRevSave=0,rightRevSave=0;
    //CASE A
    size = backward_search(csa, 0, csa.size() - 1, begin + s2, end, left, right);
    string::iterator endCaseA= begin +s2; 
    if (size != 0)
    {
        while (begin < endCaseA - 1)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = backward_search(csa, left, right, alphabet[i], leftRes, rightRes);
                if (size != 0)
                {
                    if (*(endCaseA - 1) != alphabet[i])
                        count = count + TwomissCaseASubCase(csa, leftRes, rightRes,begin,endCaseA-1,alphabet);
                     else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                    }   
                }
            }
            endCaseA--;
            left=leftSave;
            right=rightSave;
        }
    }
    //CASE B
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin, begin + s2, left, right, leftRev, rightRev);
    string::iterator beginCaseB = begin + s2;
    if (size != 0)
    {
        while (beginCaseB < end - 1)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csaReverse, leftRev, rightRev, left, right, alphabet[i], leftRevRes, rightRevRes, leftRes, rightRes);
                if (size != 0)
                {
                    if (*(beginCaseB) != alphabet[i])
                        count = count + TwomissCaseBSubCase(csa, csaReverse, leftRes,rightRes,leftRevRes,rightRevRes,beginCaseB+1,end,alphabet);
                     else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                        leftRevSave= leftRevRes;
                        rightRevSave = rightRevRes;
                    }   
                }
            }
            beginCaseB++;
            leftRev=leftRevSave;
            rightRev=rightRevSave;
            left=leftSave;
            right=rightSave;
            
        }
    }

    //CASE C
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin, begin + s1, left, right, leftRev, rightRev);
    string::iterator beginCaseC = begin+s1;
    string::iterator endCaseC = begin+s2;

    if (size != 0)
    {
        while (beginCaseC < endCaseC)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csaReverse, leftRev, rightRev, left, right, alphabet[i], leftRevRes, rightRevRes, leftRes, rightRes);
                if (size != 0)
                {
                    if (*(beginCaseC) != alphabet[i])
                    {
                        size = bidirectional_search_forward(csa, csaReverse, leftRes, rightRes, leftRevRes, rightRevRes, beginCaseC + 1, endCaseC, leftRes, rightRes, leftRevRes, rightRevRes);
                        if (size != 0)
                            count = count + TwomissCaseBSubCase (csa, csaReverse, leftRes,rightRes,leftRevRes,rightRevRes,endCaseC,end,alphabet);
                    }

                    else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                        leftRevSave= leftRevRes;
                        rightRevSave = rightRevRes;
                    }   
                }
            }
            beginCaseC++;
            leftRev=leftRevSave;
            rightRev=rightRevSave;
            left=leftSave;
            right=rightSave;
        }
    }

    //CASE D
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin + s1, begin + s2, left, right, leftRev, rightRev);
    string::iterator middle1CaseD = begin + s1;
    string::iterator middle2CaseD = begin + s2;
    if (size != 0)
    {
        while (begin < middle1CaseD)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csa , left, right, leftRev,rightRev, alphabet[i], leftRes ,rightRes, leftRevRes, rightRevRes );
                if (size != 0)
                {
                    if (*(middle1CaseD - 1) != alphabet[i])
                    {
                        size = bidirectional_search_backward(csa, csaReverse, leftRes, rightRes, leftRevRes, rightRevRes, begin, middle1CaseD - 1, leftRes, rightRes, leftRevRes, rightRevRes);
                        if (size != 0)
                            count = count + TwomissCaseBSubCase (csa, csaReverse, leftRes,rightRes,leftRevRes,rightRevRes,middle2CaseD,end,alphabet);
                    }
                    else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                        leftRevSave= leftRevRes;
                        rightRevSave = rightRevRes;
                    }   
                }
            }
            middle1CaseD--;
            leftRev=leftRevSave;
            rightRev=rightRevSave;
            left=leftSave;
            right=rightSave;
        }
    }

    return count;
}



void TwomissCaseASubCaseLocate(
    t_csa& csa,
    size_t left,
    size_t right,
    string::iterator begin,
    string::iterator end,
    const std::vector<char>& alphabet,
    int_vector<64>& results,
    size_t& index
){
    size_t leftNew,rightNew,leftSave=0,rightSave=0;
    size_t count;
    size_t leftRes=left;
    size_t rightRes=right;
    while (begin < end)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                count = backward_search(csa, leftRes, rightRes, alphabet[i],leftNew , rightNew);
                if (count != 0)
                {
                    if (*(end - 1) != alphabet[i]){
                        count = backward_search(csa, leftNew, rightNew, begin, end - 1, leftNew, rightNew);
                        for (size_t j = 0; j< count; j++)
                            results[index+j] = csa[leftNew+j];
                        index += count;
                    }
                    else{
                        leftSave=leftNew;
                        rightSave=rightNew;
                    }   
                }
            }
            end--;
            leftRes=leftSave;
            rightRes=rightSave;
        }

    
}

void TwomissCaseBSubCaseLocate(
    t_csa& csa,
    t_csa& csaReverse,
    size_t leftB,
    size_t rightB,
    size_t leftRevB,
    size_t rightRevB,
    string::iterator begin,
    string::iterator end,
    const std::vector<char>& alphabet,
    int_vector<64>& results,
    size_t& index

)
{
    size_t leftRevBRes,rightRevBRes,leftBRes,rightBRes,leftRevBSave=0,rightRevBSave=0,leftBSave=0,rightBSave=0;
    size_t count;
    size_t left=leftB;
    size_t right=rightB;
    size_t leftRev=leftRevB;
    size_t rightRev=rightRevB;
     while (begin < end)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                count = bidirectional_search(csaReverse, leftRev, rightRev, left, right, alphabet[i], leftRevBRes, rightRevBRes, leftBRes, rightBRes);
                if (count != 0)
                {
                    if (*(begin) != alphabet[i]){
                        count = bidirectional_search_forward(csa, csaReverse, leftBRes, rightBRes, leftRevBRes, rightRevBRes, begin+1, end, leftBRes, rightBRes, leftRevBRes, rightRevBRes);
                        for (size_t j = 0; j< count; j++)
                            results[index+j] = csa[leftBRes+j];
                        index += count;
                    }
                    else{
                        leftRevBSave=leftRevBRes;
                        rightRevBSave=rightRevBRes;
                        leftBSave=leftBRes;
                        rightBSave=rightBRes;
                    }    
                }
            }
            begin++;
            leftRev=leftRevBSave;
            rightRev=rightRevBSave;
            left=leftBSave;
            right=rightBSave;

        }

}

int_vector<64> locateWithTwoError( t_csa &csa,  t_csa &csaReverse, string::iterator begin,
                                                  string::iterator end,
                                                  size_t max_locate,
                                                  std::vector<char> &alphabet)
{
    size_t left;
    size_t right;
    size_t leftRev;
    size_t rightRev;
    size_t s1 = floor((double)(end-begin)/3);
    size_t s2 = end - begin - s1;
    size_t size = 0,nextIndex = 0;
    int_vector<64> results;
    results.resize(max_locate);
    size_t leftRes, rightRes, leftRevRes, rightRevRes, leftSave=0,rightSave=0,leftRevSave=0,rightRevSave=0;
    //CASE A
    size = backward_search(csa, 0, csa.size() - 1, begin + s2, end, left, right);
    string::iterator endCaseA= begin +s2; 
    if (size != 0)
    {
        while (begin < endCaseA - 1)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = backward_search(csa, left, right, alphabet[i], leftRes, rightRes);
                if (size != 0)
                {
                    if (*(endCaseA - 1) != alphabet[i])
                        TwomissCaseASubCaseLocate(csa, leftRes, rightRes,begin,endCaseA-1,alphabet,results,nextIndex);
                     else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                    }   
                }
            }
            endCaseA--;
            left=leftSave;
            right=rightSave;
        }
    }
    //CASE B
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin, begin + s2, left, right, leftRev, rightRev);
    string::iterator beginCaseB = begin + s2;
    if (size != 0)
    {
        while (beginCaseB < end - 1)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csaReverse, leftRev, rightRev, left, right, alphabet[i], leftRevRes, rightRevRes, leftRes, rightRes);
                if (size != 0)
                {
                    if (*(beginCaseB) != alphabet[i])
                        TwomissCaseBSubCaseLocate(csa, csaReverse, leftRes,rightRes,leftRevRes,rightRevRes,beginCaseB+1,end,alphabet,results,nextIndex);
                     else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                        leftRevSave= leftRevRes;
                        rightRevSave = rightRevRes;
                    }   
                }
            }
            beginCaseB++;
            leftRev=leftRevSave;
            rightRev=rightRevSave;
            left=leftSave;
            right=rightSave;
            
        }
    }

    //CASE C
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin, begin + s1, left, right, leftRev, rightRev);
    string::iterator beginCaseC = begin+s1;
    string::iterator endCaseC = begin+s2;

    if (size != 0)
    {
        while (beginCaseC < endCaseC)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csaReverse, leftRev, rightRev, left, right, alphabet[i], leftRevRes, rightRevRes, leftRes, rightRes);
                if (size != 0)
                {
                    if (*(beginCaseC) != alphabet[i])
                    {
                        size = bidirectional_search_forward(csa, csaReverse, leftRes, rightRes, leftRevRes, rightRevRes, beginCaseC + 1, endCaseC, leftRes, rightRes, leftRevRes, rightRevRes);
                        if (size != 0)
                            TwomissCaseBSubCaseLocate (csa, csaReverse, leftRes,rightRes,leftRevRes,rightRevRes,endCaseC,end,alphabet,results,nextIndex);
                    }

                    else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                        leftRevSave= leftRevRes;
                        rightRevSave = rightRevRes;
                    }   
                }
            }
            beginCaseC++;
            leftRev=leftRevSave;
            rightRev=rightRevSave;
            left=leftSave;
            right=rightSave;
        }
    }

    //CASE D
    size = bidirectional_search_forward(csa, csaReverse, 0, csa.size() - 1, 0, csaReverse.size() - 1, begin + s1, begin + s2, left, right, leftRev, rightRev);
    string::iterator middle1CaseD = begin + s1;
    string::iterator middle2CaseD = begin + s2;
    if (size != 0)
    {
        while (begin < middle1CaseD)
        {
            for (size_t i = 0; i < alphabet.size(); i++)
            {
                size = bidirectional_search(csa , left, right, leftRev,rightRev, alphabet[i], leftRes ,rightRes, leftRevRes, rightRevRes );
                if (size != 0)
                {
                    if (*(middle1CaseD - 1) != alphabet[i])
                    {
                        size = bidirectional_search_backward(csa, csaReverse, leftRes, rightRes, leftRevRes, rightRevRes, begin, middle1CaseD - 1, leftRes, rightRes, leftRevRes, rightRevRes);
                        if (size != 0)
                             TwomissCaseBSubCaseLocate (csa, csaReverse, leftRes,rightRes,leftRevRes,rightRevRes,middle2CaseD,end,alphabet,results,nextIndex);
                    }
                    else{
                        leftSave=leftRes;
                        rightSave=rightRes;
                        leftRevSave= leftRevRes;
                        rightRevSave = rightRevRes;
                    }   
                }
            }
            middle1CaseD--;
            leftRev=leftRevSave;
            rightRev=rightRevSave;
            left=leftSave;
            right=rightSave;
        }
    }

    return results;
}




int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage " << argv[0] << " text_file [max_locations] [post_context] [pre_context]" << endl;
        cout << "    This program constructs a very compact FM-index" << endl;
        cout << "    which supports count, locate, and extract queries." << endl;
        cout << "    text_file      Original text file." << endl;
        cout << "    max_locations  Maximal number of location to report." << endl;
        cout << "    post_context   Maximal length of the reported post-context." << endl;
        cout << "    pre_context    Maximal length of the pre-context." << endl;
        return 1;
    }
    size_t max_locations = 5;
    size_t post_context = 10;
    size_t pre_context = 10;
    size_t numOfErrors = 0;
    if (argc >= 3)
    {
        numOfErrors = atoi(argv[2]); //Num of errors
    }
    if (argc >= 4)
    {
        max_locations = atoi(argv[3]);
    }
    if (argc >= 5)
    {
        post_context = atoi(argv[4]);
    }
    if (argc >= 6)
    {
        pre_context = atoi(argv[5]);
    }

     int_vector<8> input;
     load_vector_from_file(input, argv[1], 1); 

     //find alphabet
     std::vector<char> alphabet;
     set<int> s( input.begin(), input.end() );
        alphabet.assign( s.begin(), s.end() );

    string index_suffix = ".fm9";
    string index_file   = string(argv[1])+index_suffix;
    t_csa fm_index;
   
    if (!load_from_file(fm_index, index_file))
    {
        ifstream in(argv[1]);
        if (!in)
        {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index " << index_file << " located. Building index now." << endl;
        construct(fm_index, argv[1], 1);    // generate index
        store_to_file(fm_index, index_file); // save it
    }

    //cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << endl;
    //cout << "Input search terms and press Ctrl-D to exit." << endl;
    //string prompt = "\e[0;32m>\e[0m ";
    //cout << prompt;
    string ReverseFile = "reversed";
    string query;
    ofstream of(ReverseFile, ofstream::binary);

    if (numOfErrors > 0)
    {
        int_vector<8> reversedInputVector(input.size());
        for (size_t i = 0; i < input.size(); i++)
            reversedInputVector[input.size()-i-1] = input[i];
        char *reversed = (char *)reversedInputVector.data();
        of.write(reversed, input.size());
        of.close();
    }

    t_csa fmRev;
    construct(fmRev, ReverseFile , 1);
    store_to_file(fm_index, "reversed"+index_suffix); // save it
    while (getline(cin, query))
    {
        size_t m = query.size();
        size_t occs = 0;
        if (numOfErrors == 0)
        {
            occs =count(fm_index, query.begin(), query.end());
        }
        if (numOfErrors == 1)
        {
            occs = countWithOneError(fm_index, fmRev, query.begin(), query.end(), alphabet);
        }
        if (numOfErrors == 2)
        {
            occs = countWithTwoError(fm_index, fmRev, query.begin(), query.end(), alphabet);
        }
        cout << query << " : " << occs << endl;
        if (occs > 0 && max_locations > 0)
        {
            int_vector<64> locations;
            cout << "Location and context of first occurrences: " << endl;
            if (numOfErrors == 0)
            {
                locations = locate(fm_index, query.begin(), query.begin() + m);
            }
            else if (numOfErrors == 1)
            {
                locations = locateWithOneError(fm_index, fmRev, query.begin(), query.begin() + m, occs, alphabet);
            }
            else if (numOfErrors == 2)
            {
                locations = locateWithTwoError(fm_index, fmRev, query.begin(), query.begin() + m, occs, alphabet);
            }
            sort(locations.begin(), locations.end());
            for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i)
            {
                cout << "  " << locations[i] << ": ";
                if (pre_extract > locations[i])
                    pre_extract = locations[i];
                
                if (locations[i] + m + post_extract > fm_index.size())
                    post_extract = fm_index.size() - locations[i] - m;
                
                auto s = extract(fm_index, locations[i] - pre_extract, locations[i] + m + post_extract - 1);
                string pre = s.substr(0, pre_extract);
                s = s.substr(pre_extract);
                if (pre.find_last_of('\n') != string::npos)
                {
                    pre = pre.substr(pre.find_last_of('\n') + 1);
                }
                cout << pre;
                //cout << "\e[1;31m";
                cout << s.substr(0, m);
                //cout << "\e[0m";
                string context = s.substr(m);
                cout << context.substr(0, context.find_first_of('\n')) << endl;
            }
        }
        //cout << prompt;
    }
    cout << endl;
}

