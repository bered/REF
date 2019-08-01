#include "Helper.h"

std::vector<std::string> splitString(const std::string &s, char delim) {
    std::stringstream ss;
    std::vector<std::string> elems;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

#ifdef _MSC_VER
#  include <intrin.h>
#  define  __builtin_popcount  __popcnt64
#endif
size_t getSupportSize(uint64_t support)
{
	return __builtin_popcount((unsigned long long)support);

	//size_t m = 0;
	//for (size_t i = 0; i < 60; i++) {
	//	if ((support & (1ull << i)) != 0)
	//		m++;
	//}
	//return m;
}

size_t getPositionOfLowestSetBit(uint64_t bitsetm) //zero based!!!!
{
    size_t m=0;
    while (true) {
        if ((bitsetm & (1ull << m)) != 0)
            break;
        m++;
    }
    return m;
}

uint64_t shiftRight(uint64_t x, size_t dimension)
{
    //see: http://blog.regehr.org/archives/1063 with n=1 and 32=dimension
    //but then the left shift has to be corrected, shift 11111... to the right for 64-dimension, then &
    return ((x>>1) | (x<<(dimension-1))) & (18446744073709551615ull>>(64-dimension));

}

uint64_t smallestRepresentation(uint64_t bitsetm, size_t dimension)
{
    uint64_t minValue=18446744073709551615ull;
    for (size_t i =0; i<dimension; i++) {
        bitsetm=shiftRight(bitsetm,dimension);
        if (bitsetm<minValue)
            minValue=bitsetm;
    }
    return minValue;
}








