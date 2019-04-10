#include "MSequenceGenerator.h"

SequenceGenerator::SequenceGenerator(unsigned int n, unsigned coefficient, unsigned seed):
    coef(coefficient),
    data(seed),
    bits(n)
{
}

int SequenceGenerator::get()
{
    auto s = 0u;
    auto t = data & coef;
    for (auto i = 0u; i < bits; ++i)
    {
        s ^= t & 1;
        t >>= 1;
    }
    data = (data >> 1) | (s << (bits - 1));
    return s;
}

void SequenceGenerator::setSeed(unsigned s)
{
    data = s;
}

void SequenceGenerator::setCoefficient(unsigned c)
{
    coef = c;
}

MSequenceGenerator::MSequenceGenerator(unsigned n):
    SequenceGenerator(n,0,0x33333333)
{
    switch (bits)
    {
    case 2:
        coef = 07;        //oct
        break;
    case 3:
        coef = 013;
        break;
    case 4:
        coef = 023;
        break;
    case 5:
        coef = 045;
        break;
    case 6:
        coef = 0103;
        break;
    case 7:
        coef = 0211;
        break;
    case 8:
        coef = 0435;
        break;
    case 9:
        coef = 01021;
        break;
    case 10:
        coef = 02011;
        break;
    case 11:
        coef = 04005;
        break;
    case 12:
        coef = 010123;
        break;
    case 13:
        coef = 020033;
        break;
    case 14:
        coef = 042103;
        break;
    case 15:
        coef = 0100003;
        break;
    case 16:
        coef = 0210013;
        break;
    case 17:
        coef = 0400011;
        break;
    case 18:
        coef = 01000201;
        break;
    case 19:
        coef = 02000047;
        break;
    case 20:
        coef = 04000011;
        break;
    case 21:
        coef = 010000005;
        break;
    case 22:
        coef = 020000003;
        break;
    case 23:
        coef = 040000041;
        break;
    case 24:
        coef = 0100000207;
        break;
    case 25:
        coef = 0200000011;
        break;
    default:
        coef = 0;
    }
    coef >>= 1;     // remove last 1.
}

InverseMSequenceGenerator::InverseMSequenceGenerator(unsigned n):
    MSequenceGenerator(n),
    square(0)
{
}

int InverseMSequenceGenerator::get()
{
    square = !square;
    return (__super::get() ^ square) * 2 - 1; // convert to -1 and 1.
}
