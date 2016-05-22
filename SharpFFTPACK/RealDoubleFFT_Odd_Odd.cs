/*
 * SharpFFTPACK is a C# version of FFTPACK (http://www.netlib.org/fftpack/).
 * SharpFFTPACK is based on jfftpack, the Java translation of fftpack (made by Baoshe Zhang).
 * It is developed as part of the SMath Studio plugin FFTPACK (made by Davide Carpi).
 * The original FFTPACK was public domain, so SharpFFTPACK is public domain too.
 * This software is in no way certified or guaranteed.
 */
using System;

namespace SharpFFTPACK
{
    /// <summary>
    /// Sine FFT transform with odd wave numbers.
    /// </summary>
    public class RealDoubleFFT_Odd_Odd : RealDoubleFFT_Even_Odd
    {
        /// <summary>
        /// Construct a wavenumber table with size n.
        /// <para>The sequences with the same size can share a wavenumber table.
        /// The prime factorization of <em>n</em> together with a tabulation of the trigonometric functions are computed and stored.</para>
        /// </summary>
        /// <param name="n">The size of a real data sequence.<br/>
        /// When <em>n</em> is a multiplication of small numbers (4, 2, 3, 5), this FFT transform is very efficient.</param>
        public RealDoubleFFT_Odd_Odd(int n)
            : base(n)
        {
        }

        /// <summary>
        /// Forward FFT transform of quarter wave data.
        /// It computes the coefficients in sine series representation with only odd wave numbers.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public new void ft(double[] x)
        {
            sinqf(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Backward FFT transform of quarter wave data.
        /// It is the unnormalized inverse transform of <em>ft</em>.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be tranformed.
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public new void bt(double[] x)
        {
            sinqb(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Forward sine FFT with odd wave numbers.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="wtable"></param>
        private void sinqf(int n, double[] x, double[] wtable)
        {
            int     k;
            double  xhold;
            int     kc, ns2;

            if (n == 1)
                return;

            ns2 = n / 2;

            for (k = 0; k < ns2; k++)
            {
                kc = n - k - 1;
                xhold = x[k];
                x[k] = x[kc];
                x[kc] = xhold;
            }

            cosqf(n, x, wtable);

            for (k = 1; k < n; k += 2)
                x[k] = -x[k];
        }

        /// <summary>
        /// Backward sine FFT with odd wave numbers.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="wtable"></param>
        private void sinqb(int n, double[] x, double[] wtable)
        {
            int     k;
            double  xhold;
            int     kc, ns2;

            if(n <= 1)
            {
                x[0] *= 4;
                return;
            }

            ns2 = n / 2;

            for (k = 1; k < n; k += 2)
                x[k] = -x[k];

            cosqb(n, x, wtable);

            for(k = 0; k < ns2; k++)
            {
                kc = n - k - 1;
                xhold = x[k];
                x[k] = x[kc];
                x[kc] = xhold;
            }
        }

//        void sinqi(int n, double[] wtable)
//        {
//            cosqi(n, wtable);
//        }
    }
}
