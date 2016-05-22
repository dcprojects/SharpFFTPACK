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
    /// Sine FFT transform of a real odd sequence.
    /// </summary>
    public class RealDoubleFFT_Odd : RealDoubleFFT_Mixed
    {
        /// <summary>
        /// <em>NormFactor</em> can be used to normalize this FFT transform.
        /// <para>This is because a call of forward transform <see cref="ft"/> followed by a call of backward transform <see cref="bt"/>
        /// will multiply the input sequence by <em>NormFactor</em></para>.
        /// </summary>
        public double NormFactor;

        private double[] WaveTable;
        private int Ndim;

        /// <summary>
        /// Construct a wavenumber table with size <em>n</em>.
        /// The sequences with the same size can share a wavenumber table.
        /// The prime factorization of <em>n</em> together with a tabulation of the trigonometric functions are computed and stored.
        /// </summary>
        /// <param name="n">The size of a real data sequence.<br/>
        /// When (<em>n</em>+1) is a multiplication of small numbers (4, 2, 3, 5), this FFT transform is very efficient.</param>
        public RealDoubleFFT_Odd(int n)
        {
            Ndim = n;
            NormFactor = 2 * (n + 1);

            int wtl = 2 * Ndim + Ndim / 2 + 3 + 15;

            if(WaveTable == null || WaveTable.Length != wtl)
                WaveTable = new double[wtl];

            sinti(Ndim, WaveTable);
        }

        /// <summary>
        /// Forward sine FFT transform.
        /// It computes the discrete sine transform of an odd sequence.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public void ft(double[] x)
        {
            sint(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Backward sine FFT transform.
        /// It is the unnormalized inverse transform of <see cref="ft"/>.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public void bt(double[] x)
        {
            sint(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Further processing of sine FFT.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="war"></param>
        /// <param name="wtable"></param>
        private void sint1(int n, double[] war, double[] wtable)
        {
            const double sqrt3 = 1.7320508075688772935274463415058723669428052538103806280558069795; // Math.Sqrt(3)

            int     modn, i, k;
            double  xhold, t1, t2;
            int     kc, np1, ns2;

            double[] wtable_p1 = new double[2 * (n + 1) + 15];

            int iw1 = n / 2;
            int iw2 = iw1 + n + 1;
            int iw3 = iw2 + n + 1;

            double[] x = new double[n + 1];

            for (i = 0; i < n; i++)
            {
                wtable[i + iw1] = war[i];
                war[i] = wtable[i + iw2];
            }

            if (n < 2)
            {
                wtable[0 + iw1] += wtable[0 + iw1];
            }
            else if (n == 2)
            {
                xhold = sqrt3 * (wtable[0 + iw1] + wtable[1 + iw1]);
                wtable[1 + iw1] = sqrt3 * (wtable[0 + iw1] - wtable[1 + iw1]);
                wtable[0 + iw1] = xhold;
            }
            else
            {
                np1 = n + 1;
                ns2 = n / 2;
                wtable[0 + iw2] = 0;

                for (k = 0; k < ns2; k++)
                {
                    kc = n - k - 1;
                    t1 = wtable[k + iw1] - wtable[kc + iw1];
                    t2 = wtable[k] * (wtable[k + iw1] + wtable[kc + iw1]);
                    wtable[k + 1 + iw2] = t1 + t2;
                    wtable[kc + 1 + iw2] = t2 - t1;
                }

                modn = n % 2;

                if (modn != 0)
                    wtable[ns2 + 1 + iw2] = 4 * wtable[ns2 + iw1];

                for (int j = 0; j < n + 1; j++)
                    wtable_p1[j] = wtable[iw1 + j];

                for (int j = n + 1; j < 2 * n + 1; j++)
                    wtable_p1[j] = war[j - (n + 1)];

                for (int j = 2 * (n + 1); j < 2 * (n + 1) + 15; j++)
                    wtable_p1[j] = wtable[iw3 + j - (2 * (n + 1))];

                for (int j = 0; j < n + 1; j++)
                    x[j] = wtable[iw2 + j];

                rfftf1(np1, x, wtable_p1, 0);

                for (int j = iw2; j < iw2 + (n + 1); j++)
                    wtable[j] = x[j - iw2];

                wtable[0 + iw1] = 0.5 * wtable[0 + iw2];

                for (i = 2; i < n; i += 2)
                {
                    wtable[i - 1 + iw1] = -wtable[i + iw2];
                    wtable[i + iw1] = wtable[i - 2 + iw1] + wtable[i - 1 + iw2];
                }

                if (modn == 0)
                    wtable[n - 1 + iw1] = -wtable[n + iw2];
            }

            for (i = 0; i < n; i++)
            {
                wtable[i + iw2] = war[i];
                war[i] = wtable[i + iw1];
            }
        }

        /// <summary>
        /// sine FFT.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <param name="wtable"></param>
        private void sint(int n, double[] x, double[] wtable)
        {
            sint1(n, x, wtable);
        }

        /// <summary>
        /// initialization of sin-FFT.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="wtable"></param>
        private void sinti(int n, double[] wtable)
        {
            int     k, ns2;
            double  dt;

            if (n <= 1)
                return;

            ns2 = n / 2;
            dt = Math.PI / (double)(n + 1);

            for(k = 0; k < ns2; k++)
                wtable[k] = 2 * Math.Sin((k + 1) * dt);

            rffti1(n + 1, wtable, ns2);
        }
    }
}
