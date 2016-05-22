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
    /// FFT transform of a complex periodic sequence.
    /// </summary>
    public class ComplexDoubleFFT : ComplexDoubleFFT_Mixed
    {
        /// <summary>
        /// <em>NormFactor</em> can be used to normalize this FFT transform.
        /// <para>This is because a call of forward transform <see cref="ft"/> followed by a call of backward transform <see cref="bt"/>
        /// will multiply the input sequence by <em>NormFactor</em>.</para>
        /// </summary>
        public double NormFactor;

        private double[] WaveTable;
        private int Ndim;
        private int Cdim;

        /// <summary>
        /// Construct a wavenumber table with size <em>n</em> for Complex FFT.
        /// <para>The sequences with the same size can share a wavenumber table.<br/>
        /// The prime factorization of <em>n</em> together with a tabulation of the trigonometric functions are computed and stored.</para>
        /// </summary>
        /// <param name="n">The size of a complex data sequence.
        /// When <i>n</i> is a multiplication of small numbers (4, 2, 3, 5), this FFT transform is very efficient.</param>
        public ComplexDoubleFFT(int n)
        {
            Ndim = n;
            Cdim = 2 * n;
            NormFactor = n;

            int wtl = 4 * Ndim + 15;

            if (WaveTable == null || WaveTable.Length != wtl)
                WaveTable = new double[wtl];

            cffti(Ndim, WaveTable);
        }

        /// <summary>
        /// Forward complex FFT transform.
        /// </summary>
        /// <param name="x">2*<em>n</em> real double data representing <em>n</em> complex double data.<br/>
        /// As an input parameter, <em>x</em> is an array of 2*<em>n</em> real data representing <em>n</em> complex data.<br/>
        /// As an output parameter, <em>x</em> represents <em>n</em> FFT'd complex data.<br/>
        /// Their relation as follows:<br/>
        /// x[2*i] is the real part of <em>i</em>-th complex data;<br/>
        /// x[2*i+1] is the imaginary part of <em>i</em>-th complex data.
        /// </param>
        public void ft(double[] x)
        {
            if (x.Length != Cdim)
                throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");

            cfftf(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Forward complex FFT transform.
        /// </summary>
        /// <param name="x">An array of <em>n</em> Complex data.</param>
        public void ft(Complex[] x)
        {
            if (x.Length != Ndim)
                throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");

            double[] y = new double[Cdim];

            for (int i = 0; i < Ndim; i++)
            {
                y[2 * i] = x[i].Re;
                y[2 * i + 1] = x[i].Im;
            }

            cfftf(Ndim, y, WaveTable);

            for (int i = 0; i < Ndim; i++)
            {
                x[i].Re = y[2 * i];
                x[i].Im = y[2 * i + 1];
            }
        }

        /// <summary>
        /// Backward complex FFT transform. It is the unnormalized inverse transform of <see cref="ft(double[])"/>.
        /// </summary>
        /// <param name="x">2*<em>n</em> real double data representing <em>n</em> complex double data.
        /// As an input parameter, <em>x</em> is an array of 2*<em>n</em> real data representing <em>n</em> complex data.<br/>
        /// As an output parameter, <em>x</em> represents <em>n</em> FFT'd complex data.<br/>
        /// Their relation as follows:<br/>
        /// x[2*<em>i</em>] is the real part of <em>i</em>-th complex data;<br/>
        /// x[2*<em>i</em>+1] is the imaginary part of <em>i</em>-th complex data.
        /// </param>
        public void bt(double[] x)
        {
            if (x.Length != Cdim)
                throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");

            cfftb(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Backward complex FFT transform. It is the unnormalized inverse transform of <see cref="ft(Complex1D)"/>.
        /// </summary>
        /// <param name="x">an array of <em>n</em> Complex data.</param>
        public void bt(Complex[] x)
        {
            if (x.Length != Ndim)
                throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");

            double[] y = new double[Cdim];

            for (int i = 0; i < Ndim; i++)
            {
                y[2 * i] = x[i].Re;
                y[2 * i + 1] = x[i].Im;
            }

            cfftb(Ndim, y, WaveTable);

            for (int i = 0; i < Ndim; i++)
            {
                x[i].Re = y[2 * i];
                x[i].Im = y[2 * i + 1];
            }
        }
    }
}
