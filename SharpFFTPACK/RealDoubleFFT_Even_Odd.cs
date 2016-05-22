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
    /// Cosine FFT transform with odd wave numbers.
    /// </summary>
    public class RealDoubleFFT_Even_Odd : RealDoubleFFT_Mixed
    {
        /// <summary>
        /// <em>NormFactor</em> can be used to normalize this FFT transform.
        /// <para>This is because a call of forward transform <see cref="ft"/> followed by a call of backward transform <see cref="bt"/>
        /// will multiply the input sequence by <em>NormFactor</em></para>.
        /// </summary>
        public double NormFactor;

        protected double[] WaveTable;
        protected int Ndim;

        /// <summary>
        /// Construct a wavenumber table with size <em>n</em>.
        /// The sequences with the same size can share a wavenumber table.
        /// The prime factorization of <em>n</em> together with a tabulation of the trigonometric functions are computed and stored.
        /// </summary>
        /// <param name="n">The size of a real data sequence.<br/>
        /// When <em>n</em> is a multiplication of small numbers (4, 2, 3, 5), this FFT transform is very efficient.</param>
        public RealDoubleFFT_Even_Odd(int n)
        {
            Ndim = n;
            NormFactor = 4 * n;

            int wtl = 3 * Ndim + 15;

            if(WaveTable == null || WaveTable.Length != wtl)
                WaveTable = new double[wtl];

            cosqi(Ndim, WaveTable);
        }

        /// <summary>
        /// Forward FFT transform of quarter wave data.
        /// It computes the coefficients in cosine series representation with only odd wave numbers.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public void ft(double[] x)
        {
            cosqf(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Backward FFT transform of quarter wave data.
        /// It is the unnormalized inverse transform of <em>ft</em>.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be tranformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public void bt(double[] x)
        {
            cosqb(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Further processing of forward cos-FFT with odd wave numbers.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="wtable"></param>
        private void cosqf1(int n, double[] x, double[] wtable)
        {
            int     modn, i, k;
            int     kc, np2, ns2;
            double  xim1;

            ns2 = (n + 1) / 2;
            np2 = n + 2;

            for (k = 1; k < ns2; k++)
            {
                kc = n - k;
                wtable[k + n] = x[k] + x[kc];
                wtable[kc + n] = x[k] - x[kc];
            }

            modn = n % 2;

            if (modn == 0)
                wtable[ns2 + n] = x[ns2] + x[ns2];

            for (k = 1; k < ns2; k++)
            {
                kc = n - k;
                x[k] = wtable[k - 1] * wtable[kc + n] + wtable[kc - 1] * wtable[k + n];
                x[kc] = wtable[k - 1] * wtable[k + n] - wtable[kc - 1] * wtable[kc + n];
            }

            if(modn == 0)
                x[ns2] = wtable[ns2 - 1] * wtable[ns2 + n];

            rfftf1(n, x, wtable, n);

            for (i = 2; i < n; i += 2)
            {
                xim1 = x[i - 1] - x[i];
                x[i] = x[i - 1] + x[i];
                x[i - 1] = xim1;
            }
        }

        /// <summary>
        /// Further processing of backward cos-FFT with odd wave numbers.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="wtable"></param>
        private void cosqb1(int n, double[] x, double[] wtable)
        {
            int     modn, i, k;
            int     kc, ns2;
            double  xim1;

            ns2 = (n + 1) / 2;

            for (i = 2; i < n; i += 2)
            {
                xim1 = x[i - 1] + x[i];
                x[i] -= x[i - 1];
                x[i - 1] = xim1;
            }

            x[0] += x[0];
            modn = n % 2;

            if (modn == 0)
                x[n - 1] += x[n - 1];

            rfftb1(n, x, wtable, n);

            for (k = 1; k < ns2; k++)
            {
                kc = n - k;
                wtable[k + n] = wtable[k - 1] * x[kc] + wtable[kc - 1] * x[k];
                wtable[kc + n] = wtable[k - 1] * x[k] - wtable[kc - 1] * x[kc];
            }

            if (modn == 0)
                x[ns2] = wtable[ns2 - 1] * (x[ns2] + x[ns2]);

            for (k = 1; k < ns2; k++)
            {
                kc = n - k;
                x[k] = wtable[k + n] + wtable[kc + n];
                x[kc] = wtable[k + n] - wtable[kc + n];
            }

            x[0] += x[0];
        }

        /// <summary>
        /// Forward cosine FFT with odd wave numbers.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="wtable"></param>
        protected void cosqf(int n, double[] x, double[] wtable)
        {
            const double sqrt2 = 1.4142135623730950488016887242096980785696718753769480731766797380; // Math.Sqrt(2)

            double  tsqx;

            if (n < 2)
                return;

            if (n == 2)
            {
                tsqx = sqrt2 * x[1];
                x[1] = x[0] - tsqx;
                x[0] += tsqx;
            }
            else
            {
                cosqf1(n, x, wtable);
            }
        }

        /// <summary>
        /// Backward cosine FFT with odd wave numbers.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="wtable"></param>
        protected void cosqb(int n, double[] x, double[] wtable)
        {
            const double tsqrt2 = 2.8284271247461900976033774484193961571393437507538961463533594760; // 2 * Math.Sqrt(2)

            double  x1;

            if (n < 2)
            {
                x[0] *= 4;
            }
            else if (n == 2)
            {
                x1 = 4 * (x[0] + x[1]);
                x[1] = tsqrt2 * (x[0] - x[1]);
                x[0] = x1;
            }
            else
            {
                cosqb1(n, x, wtable);
            }
        }

        /// <summary>
        /// Initialization of cosine FFT with odd wave numbers.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="wtable"></param>
        private void cosqi(int n, double[] wtable)
        {
            const double pih = Math.PI / 2D;

            int     k;
            double  dt;

            dt = pih / (double)n;

            for (k = 0; k < n; k++)
                wtable[k] = Math.Cos((k + 1) * dt);

            rffti1(n, wtable, n);
        }
    }
}
