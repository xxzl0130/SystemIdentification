#pragma once

/**
 * \brief Simple 01 sequence generator.
 */
class SequenceGenerator
{
public:
    SequenceGenerator(unsigned int n, unsigned int coefficient = 07,unsigned int seed = 1);
    virtual ~SequenceGenerator() = default;

    // Get next bit of sequence.
    virtual int get();
    // Set the seed (initial state) of sequence.
    virtual void setSeed(unsigned int s);
    // Set the coefficient of polynomial.
    virtual void setCoefficient(unsigned int c);
protected:
    unsigned int coef;      // Coefficient of Primitive polynomial.
    unsigned int data;      // Sequence data.
    unsigned int bits;      // Cycle bits of sequence.
};

/**
 * \brief M sequence generator, using primitive polynomial to generator maximum cycle of 01 sequence.
 *Sequence cycle will be (2^n-1)bits. Currently support up to n = 25.
 */
class MSequenceGenerator : public SequenceGenerator
{
public:
    MSequenceGenerator(unsigned int n = 2);
    virtual ~MSequenceGenerator() = default;

protected:
};

/**
 * \brief Inverse M sequence generator.
 */
class InverseMSequenceGenerator : public MSequenceGenerator
{
public:
    InverseMSequenceGenerator(unsigned int n = 2);
    virtual ~InverseMSequenceGenerator() = default;

    // Get next bit of inverse m sequence
    virtual int get() override;
protected:
    unsigned int square;    //Square wave data.
};