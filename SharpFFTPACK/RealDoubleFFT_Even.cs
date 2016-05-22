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
    /// Cosine FFT transform of a real even sequence.
    /// </summary>
    public class RealDoubleFFT_Even : RealDoubleFFT_Mixed
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
        /// <param name="n">The size of a real data sequence.
        /// When (<em>n</em>-1) is a multiplication of small numbers (4, 2, 3, 5), this FFT transform is very efficient.</param>
        public RealDoubleFFT_Even(int n)
        {
            Ndim = n;
            NormFactor = 2 * (n - 1);

            int wtl = 3 * Ndim + 15;

            if (WaveTable == null || WaveTable.Length != wtl)
                WaveTable = new double[wtl];

            costi(Ndim, WaveTable);
        }

        /// <summary>
        /// Forward cosine FFT transform.
        /// It computes the discrete sine transform of an odd sequence.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public void ft(double[] x)
        {
            cost(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Backward cosine FFT transform.
        /// It is the unnormalized inverse transform of <see cref="ft"/>.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        public void bt(double[] x)
        {
            cost(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Cosine FFT. Backward and forward cos-FFT are the same.
        /// </summary>
        private void cost(int n, double[] x, double[] wtable)
        {
            int     modn, i, k;
            double  c1, t1, t2;
            int     kc;
            double  xi;
            int     nm1;
            double  x1h;
            int     ns2;
            double  tx2, x1p3, xim2;

            nm1 = n - 1;
            ns2 = n / 2;

            if (n - 2 < 0)
                return;

            if (n == 2)
            {
                x1h = x[0] + x[1];
                x[1] = x[0] - x[1];
                x[0] = x1h;
            }
            else if (n == 3)
            {
                x1p3 = x[0] + x[2];
                tx2 = x[1] + x[1];
                x[1] = x[0] - x[2];
                x[0] = x1p3 + tx2;
                x[2] = x1p3 - tx2;
            }
            else
            {
                c1 = x[0] - x[n - 1];
                x[0] += x[n - 1];

                for (k = 1; k < ns2; k++)
                {
                    kc = nm1 - k;
                    t1 = x[k] + x[kc];
                    t2 = x[k] - x[kc];
                    c1 += wtable[kc] * t2;
                    t2 = wtable[k] * t2;
                    x[k] = t1 - t2;
                    x[kc] = t1 + t2;
                }

                modn = n % 2;

                if (modn !=0)
                    x[ns2] += x[ns2];

                rfftf1(nm1, x, wtable, n);
                xim2 = x[1];
                x[1] = c1;

                for (i = 3; i < n; i += 2)
                {
                    xi = x[i];
                    x[i] = x[i - 2] - x[i - 1];
                    x[i - 1] = xim2;
                    xim2 = xi;
                }

                if (modn != 0)
                    x[n - 1] = xim2;
            }
        }

        /// <summary>
        /// initialization of cos-FFT.
        /// </summary>
        /// <param name="n"></param>
        /// <param name="wtable"></param>
        private void costi(int n, double[] wtable)
        {
            int     k, kc, ns2;
            double  dt;

            if (n <= 3)
                return;

            ns2 = n / 2;
            dt = Math.PI / (double)(n - 1);

            for (k = 1; k < ns2; k++)
            {
                kc = n - k - 1;
                wtable[k] = 2 * Math.Sin(k * dt);
                wtable[kc] = 2 * Math.Cos(k * dt);
            }

            rffti1(n - 1, wtable, n);
        }
    }
}
