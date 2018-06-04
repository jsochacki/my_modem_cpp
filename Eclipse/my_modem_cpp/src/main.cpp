
#include <iostream>
#include <random>
#include <complex>

//can be set to float , double, or long double only
using simulation_precision = double;

//Equivalent to
//std::mt19937_64 rng;
std::mersenne_twister_engine<uint_fast64_t,
                             64, 312, 156, 31,
                             0xb5026f5aa96619e9ULL, 29,
                             0x5555555555555555ULL, 17,
                             0x71d67fffeda60000ULL, 37,
                             0xfff7eee000000000ULL, 43,
                             6364136223846793005ULL> rng;

std::uniform_real_distribution<simulation_precision> mt_uniform(0.0, 1.0);

template <typename Type>
std::complex<Type> complex_gaussian_point_noise(void);

template <typename Type>
inline std::complex<Type> inline_complex_gaussian_point_noise(void);

template <typename Type>
Type standard_deviation(unsigned int, Type *);

template <typename Type>
std::complex<Type> standard_deviation(unsigned int, std::complex<Type> *);

template <typename Type>
Type mean_squared_signal_power(unsigned int, std::complex<Type> *);

template <typename Type>
Type mean_signal_power(unsigned int , std::complex<Type> *);

template <typename Type>
void AWGN_Generator(unsigned int,
                    std::complex<Type> *,
                    Type,
                    unsigned int);

template <typename Type>
void AGC(unsigned int, std::complex<Type> *, Type);

int main(int argc, char **argv)
{


   //Not what it looks like
   //std::seed_seq sseq{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
   //rng.seed(sseq);

   rng.seed(0);

   unsigned int nn = 1000000;
   std::complex<simulation_precision> *values = new std::complex<simulation_precision>[nn]();
   simulation_precision *dvalues = new simulation_precision[nn]();
   for(int n = 0; n < nn; n++)
   {
      values[n] = complex_gaussian_point_noise<simulation_precision>();
   }
   for(int n = 0; n < nn; n += 2)
   {
      dvalues[n] = complex_gaussian_point_noise<simulation_precision>().real();
      dvalues[n+1] = complex_gaussian_point_noise<simulation_precision>().imag();
   }

   std::cout << standard_deviation(nn, values) << std::endl;
   std::cout << standard_deviation(nn, dvalues) << std::endl;

   std::cout << mean_squared_signal_power(nn, values) << std::endl;
   std::cout << mean_signal_power(nn, values) << std::endl;

   AWGN_Generator<simulation_precision>(nn, values, 3.0, 2);
   AGC<simulation_precision>(nn, values, 1.0);
   std::cout << standard_deviation(nn, values) << std::endl;
   std::cout << mean_squared_signal_power(nn, values) << std::endl;
   std::cout << mean_signal_power(nn, values) << std::endl;

   delete[] values;
   delete[] dvalues;
}

//supports float, double, or long double
template <typename Type>
std::complex<Type> complex_gaussian_point_noise(void)
{
   std::complex<Type> urv; Type r, n;

   do
   {
      urv = std::complex<Type>( (mt_uniform(rng) * 2) - 1, (mt_uniform(rng) * 2) - 1);
      r = std::pow(std::abs(urv), 2);
   }
   while ( ( r >= 1) || (r == 0 ) );

   n = std::sqrt(-2 * std::log(r) / r);
   return std::complex<Type>(urv * n);
}

//supports float, double, or long double
template <typename Type>
inline std::complex<Type> inline_complex_gaussian_point_noise(void)
{
   std::complex<Type> urv; Type r, n;

   do
   {
      urv = std::complex<Type>( (mt_uniform(rng) * 2) - 1, (mt_uniform(rng) * 2) - 1);
      r = std::pow(std::abs(urv), 2);
   }
   while ( ( r >= 1) || (r == 0 ) );

   n = std::sqrt(-2 * std::log(r) / r);
   return std::complex<Type>(urv * n);
}

//supports float, double, or long double
template <typename Type>
std::complex<Type> standard_deviation(unsigned int input_length,
                                        std::complex<Type> *input_values)
{
   unsigned int n = 0;
   std::complex<Type> mean(0.0, 0.0), stdev(0.0, 0.0), tvar(0.0, 0.0);

   for(n = 0; n < input_length; n++)
   {
      mean += input_values[n];
   }

   mean.real(mean.real() / input_length);
   mean.imag(mean.imag() / input_length);

   for(n = 0; n < input_length; n++)
   {
      tvar = input_values[n] - mean;
      tvar.real(std::pow(tvar.real(), 2));
      tvar.imag(std::pow(tvar.imag(), 2));
      stdev += tvar;
   }
   stdev.real(stdev.real() / input_length);
   stdev.imag(stdev.imag() / input_length);

   return stdev;
}

//supports anything std::pow does
template <typename Type>
Type standard_deviation(unsigned int input_length, Type *input_values)
{
   unsigned int n = 0;
   Type mean(0), stdev(0);

   for(n = 0; n < input_length; n++)
   {
      mean += input_values[n];
   }

   mean /= input_length;

   for(n = 0; n < input_length; n++)
   {
      stdev += std::pow(input_values[n] - mean, 2);
   }
   stdev /= input_length;

   return stdev;
}

//supports float, double, or long double
template <typename Type>
Type mean_squared_signal_power(unsigned int signal_length, std::complex<Type> *symbol_stream)
{
    unsigned int n;
    Type sss;

    sss = 0, n = 0;

    for(n = 0; n < signal_length; n++)
    {
        sss += std::pow( std::abs(symbol_stream[n]) , 2);
    }

    return ((1 / static_cast<Type>(signal_length)) * sss);
}

//supports float, double, or long double
template <typename Type>
Type mean_signal_power(unsigned int signal_length, std::complex<Type> *symbol_stream)
{
    unsigned int n;
    std::complex<Type> sss(0.0, 0.0);

    n = 0;

    for(n = 0; n < signal_length; n++)
    {
        sss += symbol_stream[n];
    }

    return ((1 / static_cast<Type>(signal_length)) * std::abs(sss));
}

//This function does change the noise level based on the signal level so that
//It is correct for any signal power
template <typename Type>
void AWGN_Generator(unsigned int signal_length, std::complex<Type> *symbol_stream, Type EbN0_dB, unsigned int bits_per_symbol)
{
   Type EbN0, SNR, mssp, N0, stdev;

   EbN0 = std::pow(10, EbN0_dB / 10);
   SNR = EbN0 * bits_per_symbol;
   mssp = mean_squared_signal_power<Type>(signal_length, symbol_stream);
   N0 = mssp / SNR;
   stdev = std::sqrt(N0 / 2);

   //Generate Gaussian Noise with appropriate variance
   unsigned int n = 0;
   std::complex<Type> *noise;
   std::complex<Type> stdev_value(0.0, 0.0);
   std::complex<Type> scaling_factor(0.0, 0.0);

   noise = new std::complex<Type>[signal_length]();

   for(n = 0; n < signal_length; n++)
   {
      noise[n] = inline_complex_gaussian_point_noise<Type>();
   }

   stdev_value = standard_deviation( signal_length, noise);

   scaling_factor.real( std::sqrt( std::pow(stdev, 2) / std::pow(stdev_value.real(), 2)));
   scaling_factor.imag( std::sqrt( std::pow(stdev, 2) / std::pow(stdev_value.imag(), 2)));

   // Add gaussian noise and scale at the same time for speed

   for (n = 0; n < signal_length; n++)
   {
     symbol_stream[n].real( symbol_stream[n].real() +
                             (noise[n].real() * scaling_factor.real()));
     symbol_stream[n].imag( symbol_stream[n].imag() +
                             (noise[n].imag() * scaling_factor.imag()));
   }

   delete[] noise;
}

template <typename Type>
void AGC(unsigned int signal_length, std::complex<Type> *symbol_stream, Type desired_sum_squared_power)
{
    unsigned int n;
    Type A;

    A = std::sqrt( mean_squared_signal_power<Type>( signal_length, symbol_stream) / desired_sum_squared_power);

    for(n = 0; n < signal_length; n++)
    {
        symbol_stream[n] /= A;
    }
}

//TODO Finish up converting these
//void modulator(int modulation, int qam_length, COMPLEX *qam_symbols,
//
//               bool * encoded_bits, COMPLEX *channel_qam, bool rc_interleave)
//{
//
//
//    int i, b, m, s;
//
//    /* Map coded bits to constellation symbols */
//
//    if (rc_interleave)
//    {
//        for (b = 0;  b < qam_length;  b++)
//        {
//            for (s = 0, m = 0;  m < modulation; m++)
//            {
//               int m2 = m;
//               s = s + s + encoded_bits[b+m2*qam_length];
//            }
//        channel_qam[b] = qam_symbols[s];
//        }
//    }
//    else
//    {
//        for (b = 0;  b < qam_length;  b++)
//        {
//            for (s = 0, m = 0;  m < modulation; m++)
//            {
//               int m2 = m;
//               s = s + s + encoded_bits[b*modulation+m2];
//            }
//        channel_qam[b] = qam_symbols[s];
//        }
//    }
//
//
//}
//
///*************************************************************************************/
//void demodulator (int modulation, int qam_length,
//
//                 QUANT *LLR, COMPLEX *channel_qam, QUANT *channel_soft,
//
//                 int iq_limit, int siv_limit, double iq_quant_ratio, double siv_quant_ratio, bool rc_interleave)
//{
//
//    QUANT i, q;
//    int b, m, j;
//
//    for (b = 0;  b < qam_length;  b++)
//    {
//        i = quantize(channel_qam[b].i, iq_quant_ratio, iq_limit);
//        q = quantize(channel_qam[b].q, iq_quant_ratio, iq_limit);
//
//        for (m=0; m<modulation; m++)
//        {
//            if (rc_interleave)
//            {
//               int m2 = modulation-1-m;
//               channel_soft[b+m*qam_length] = LLR[(i+iq_limit)*2*iq_limit*modulation+(q+iq_limit)*modulation+m2];
//            }
//            else
//            {
//               int m2 = modulation-1-m;
//               channel_soft[b*modulation+m] = LLR[(i+iq_limit)*2*iq_limit*modulation+(q+iq_limit)*modulation+m2];
//            }
//        }
//    }
//
//}
//
//
///*************************************************************************************/
//void qpsk_maximum_liklihood_hard_decision_decoder (int modulation, int qam_length, COMPLEX *qam_symbols,
//
//                                                   DIGITAL *encoded_bits, COMPLEX *channel_qam, QUANT * channel_hard)
//{
//
//    //qam_symbols is the complex value for the different symbols
//    //so you will use qam_symbols.i for the i value and qam_symbols.q for the q value
//    //channel_qam is the complex value of the symbol after the channel impairments are added
//    //qam_symbols are the complex symbols
//    //qam_symbols[0] corresponds to binary MSB 0 0 LSB
//    //qam_symbols[1] corresponds to binary MSB 0 1 LSB
//    //qam_symbols[2] corresponds to binary MSB 1 0 LSB
//    //qam_symbols[3] corresponds to binary MSB 1 1 LSB
//
////    FILE *fname, *fname2;
////
////    fname = fopen("qpskdata.txt", "w");
////    fname2 = fopen("encodedbits.txt", "w");
//
//    COMPLEX temp;
//    double current_magnitude, smallest_magnitude;
//    int closest;
//    int s, d, m, mtemp;
//
//    //may be slightly biased as it pays no special treatment to symbols of equal
//    //magnitude and will always be equal to the lower bit value when there
//    //are hard decisions that are a tie
//    for (s = 0;  s < qam_length;  s++)
//    {
////        fprintf(fname, "%lf , %lf\n", channel_qam[s].i, channel_qam[s].q);
////        fprintf(fname2, "%d\n%d\n", encoded_bits[s*2], encoded_bits[(s*2) + 1]);
//        temp.i = 0, temp.q = 0;
//        current_magnitude = 0, smallest_magnitude = 100;
//        closest = 0;
//
//        for (d = 0; d < (1 << modulation); d++)
//        {
//            temp.i = channel_qam[s].i - qam_symbols[d].i;
//            temp.q = channel_qam[s].q - qam_symbols[d].q;
//            current_magnitude = sqrt(pow(temp.i, 2) + pow(temp.q, 2));
//            if (current_magnitude < smallest_magnitude)
//            {
//                closest = d;
//                smallest_magnitude = current_magnitude;
//            }
//        }
//
//        for (m = 0; m < modulation; m++)
//        {
//               mtemp = modulation-1-m;
//               channel_hard[s*modulation+mtemp] = (QUANT) ((closest & (m + 1)) > 0);
//        }
//    }
////    fclose(fname);
////    fclose(fname2);
//
//}





/*


//C only implementations

typedef double ANALOG;

typedef struct  COMPLEX {

   ANALOG i;

   ANALOG q;

};


void modulator(int modulation, int qam_length, COMPLEX *qam_symbols,

               bool * encoded_bits, COMPLEX *channel_qam, bool rc_interleave)
{


    int i, b, m, s;

     Map coded bits to constellation symbols

    if (rc_interleave)
    {
        for (b = 0;  b < qam_length;  b++)
        {
            for (s = 0, m = 0;  m < modulation; m++)
            {
               int m2 = m;
               s = s + s + encoded_bits[b+m2*qam_length];
            }
        channel_qam[b] = qam_symbols[s];
        }
    }
    else
    {
        for (b = 0;  b < qam_length;  b++)
        {
            for (s = 0, m = 0;  m < modulation; m++)
            {
               int m2 = m;
               s = s + s + encoded_bits[b*modulation+m2];
            }
        channel_qam[b] = qam_symbols[s];
        }
    }


}


void demodulator (int modulation, int qam_length,

                 QUANT *LLR, COMPLEX *channel_qam, QUANT *channel_soft,

                 int iq_limit, int siv_limit, double iq_quant_ratio, double siv_quant_ratio, bool rc_interleave)
{

    QUANT i, q;
    int b, m, j;

    for (b = 0;  b < qam_length;  b++)
    {
        i = quantize(channel_qam[b].i, iq_quant_ratio, iq_limit);
        q = quantize(channel_qam[b].q, iq_quant_ratio, iq_limit);

        for (m=0; m<modulation; m++)
        {
            if (rc_interleave)
            {
               int m2 = modulation-1-m;
               channel_soft[b+m*qam_length] = LLR[(i+iq_limit)*2*iq_limit*modulation+(q+iq_limit)*modulation+m2];
            }
            else
            {
               int m2 = modulation-1-m;
               channel_soft[b*modulation+m] = LLR[(i+iq_limit)*2*iq_limit*modulation+(q+iq_limit)*modulation+m2];
            }
        }
    }

}



void qpsk_maximum_liklihood_hard_decision_decoder (int modulation, int qam_length, COMPLEX *qam_symbols,

                                                   DIGITAL *encoded_bits, COMPLEX *channel_qam, QUANT * channel_hard)
{

    //qam_symbols is the complex value for the different symbols
    //so you will use qam_symbols.i for the i value and qam_symbols.q for the q value
    //channel_qam is the complex value of the symbol after the channel impairments are added
    //qam_symbols are the complex symbols
    //qam_symbols[0] corresponds to binary MSB 0 0 LSB
    //qam_symbols[1] corresponds to binary MSB 0 1 LSB
    //qam_symbols[2] corresponds to binary MSB 1 0 LSB
    //qam_symbols[3] corresponds to binary MSB 1 1 LSB

//    FILE *fname, *fname2;
//
//    fname = fopen("qpskdata.txt", "w");
//    fname2 = fopen("encodedbits.txt", "w");

    COMPLEX temp;
    double current_magnitude, smallest_magnitude;
    int closest;
    int s, d, m, mtemp;

    //may be slightly biased as it pays no special treatment to symbols of equal
    //magnitude and will always be equal to the lower bit value when there
    //are hard decisions that are a tie
    for (s = 0;  s < qam_length;  s++)
    {
//        fprintf(fname, "%lf , %lf\n", channel_qam[s].i, channel_qam[s].q);
//        fprintf(fname2, "%d\n%d\n", encoded_bits[s*2], encoded_bits[(s*2) + 1]);
        temp.i = 0, temp.q = 0;
        current_magnitude = 0, smallest_magnitude = 100;
        closest = 0;

        for (d = 0; d < (1 << modulation); d++)
        {
            temp.i = channel_qam[s].i - qam_symbols[d].i;
            temp.q = channel_qam[s].q - qam_symbols[d].q;
            current_magnitude = sqrt(pow(temp.i, 2) + pow(temp.q, 2));
            if (current_magnitude < smallest_magnitude)
            {
                closest = d;
                smallest_magnitude = current_magnitude;
            }
        }

        for (m = 0; m < modulation; m++)
        {
               mtemp = modulation-1-m;
               channel_hard[s*modulation+mtemp] = (QUANT) ((closest & (m + 1)) > 0);
        }
    }
//    fclose(fname);
//    fclose(fname2);

}


void AGC(int qam_length, COMPLEX *qam_symbols, double desired_sum_squared_power)
{
    int n;
    double A;

    A = sqrt( mean_squared_signal_power( qam_length, qam_symbols) / desired_sum_squared_power);

    for(n = 0; n < qam_length; n++)
    {
        qam_symbols[n].i = qam_symbols[n].i / A;
        qam_symbols[n].q = qam_symbols[n].q / A;
    }
}


//This function does change the noise level based on the signal level so that
//It is correct for any signal power
void AWGN_Generator(int qam_length, COMPLEX *qam_symbols, double EbN0_dB, int modulation)
{

//    FILE *noisefile;
//
//    noisefile = fopen("noisefile.txt", "w");

    double EbN0, SNR, mssp, N0, stdev;

    EbN0 = pow(10, EbN0_dB / 10);
    SNR = EbN0 * modulation;
    mssp = mean_squared_signal_power( qam_length, qam_symbols);
    N0 = mssp / SNR;
    stdev = sqrt(N0 / 2);

    //Generate Gaussian Noise with appropriate variance
    int n = 0;
    double r, s;
    COMPLEX *noise;
    COMPLEX stdev_value;
    COMPLEX scaling_factor;

    stdev_value.i = 0, stdev_value.q = 0;
    scaling_factor.i = 0, scaling_factor.q = 0;

    noise = (COMPLEX*) calloc(qam_length, sizeof(COMPLEX));

    for(n = 0; n < qam_length; n++)
    {
        do
        {
            noise[n].i = ((mt_uniform()*2.0 ) - 1);
            noise[n].q = ((mt_uniform()*2.0 ) - 1);
            r = (noise[n].i * noise[n].i) + (noise[n].q * noise[n].q);
        }
        while ((r >= 1) || (r == 0));
        s = sqrt(-2 * log(r) / r);
        noise[n].i = noise[n].i * s;
        noise[n].q = noise[n].q * s;
    }

    stdev_value = complex_standard_deviation( qam_length, noise);

    scaling_factor.i = sqrt(pow(stdev, 2) / pow(stdev_value.i, 2));
    scaling_factor.q = sqrt(pow(stdev, 2) / pow(stdev_value.q, 2));

//    for(n = 0; n < qam_length; n++)
//    {
//        noise[n].i *= scaling_factor.i;
//        noise[n].q *= scaling_factor.q;
//    }

    // Add gaussian noise and scale at the same time for speed

    for (n = 0; n < qam_length; n++) {
        qam_symbols[n].i += (noise[n].i * scaling_factor.i);
        qam_symbols[n].q += (noise[n].q * scaling_factor.q);
//        fprintf(noisefile,"%lf,%lf\n",(noise[n].i * scaling_factor.i),(noise[n].q * scaling_factor.q));
    }
//    fclose(noisefile);
    free(noise);
}


double mean_squared_signal_power(int qam_length, COMPLEX *qam_symbols)
{
    int n;
    double sss;

    sss = 0, n = 0;

    for(n = 0; n < qam_length; n++)
    {
        sss += pow(qam_symbols[n].i, 2) + pow(qam_symbols[n].q, 2);
    }

    return ((1 / (double)qam_length) * sss);
}


double double_standard_deviation(int input_length, double *input_values)
{
    int n = 0;
    double mean = 0;
    double stdev = 0;

    for(n = 0; n < input_length; n++)
    {
        mean += input_values[n];
    }
    mean = mean / input_length;

    for(n = 0; n < input_length; n++)
    {
        stdev += pow(input_values[n] - mean, 2);
    }
    stdev = stdev / input_length;

    return sqrt(stdev);

}


COMPLEX complex_standard_deviation(int input_length, COMPLEX *input_values)
{
    int n = 0;
    COMPLEX mean;
    COMPLEX stdev;

    mean.i = 0, mean.q = 0;
    stdev.i = 0, stdev.q = 0;

    for(n = 0; n < input_length; n++)
    {
        mean.i += input_values[n].i;
        mean.q += input_values[n].q;
    }
    mean.i = mean.i / input_length;
    mean.q = mean.q / input_length;

    for(n = 0; n < input_length; n++)
    {
        stdev.i += pow(input_values[n].i - mean.i, 2);
        stdev.q += pow(input_values[n].q - mean.q, 2);
    }
    stdev.i = sqrt(stdev.i / input_length);
    stdev.q = sqrt(stdev.q / input_length);

    return stdev;

}


//Has a standard deviation of 1 and a mean of 0
COMPLEX complex_gaussian_point_noise(void)
{
    COMPLEX urv, rv;
    double r, n;

    do
    {
      urv.i = (mt_uniform()*2) - 1;
      urv.q = (mt_uniform()*2) - 1;
      r = (urv.i * urv.i)+(urv.q * urv.q);
    }
    while ((r >= 1) || (r == 0));

    n = sqrt(-2 * log(r) / r);
    rv.i = urv.i * n;
    rv.q = urv.q * n;

    return (rv);
}

*/

