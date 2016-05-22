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
    /// FFT transform of a real periodic sequence.
    /// </summary>
    public class RealDoubleFFT : RealDoubleFFT_Mixed
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
        /// <para>The sequences with the same size can share a wavenumber table.<br/>
        /// The prime factorization of <em>n</em> together with a tabulation of the trigonometric functions are computed and stored.</para>
        /// </summary>
        /// <param name="n">The size of a real data sequence.
        /// When <em>n</em> is a multiplication of small numbers (4, 2, 3, 5), this FFT transform is very efficient.</param>
        public RealDoubleFFT(int n)
        {
            Ndim = n;
            NormFactor = n;

            int wtl = 2 * Ndim + 15;

            if (WaveTable == null || WaveTable.Length != wtl)
                WaveTable = new double[wtl];

            rffti(Ndim, WaveTable);
        }

        /// <summary>
        /// Forward real FFT transform. It computes the discrete transform of a real data sequence.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// <para>After FFT, <em>x</em> contains the transform coefficients used to construct <em>n</em> complex FFT coefficients.
        /// The real part of the first complex FFT coefficients is <em>x</em>[0]; its imaginary part is 0.</para>
        /// <para>If <em>n</em> is even set <em>m</em> = <em>n</em>/2, if <em>n</em> is odd set <em>m</em> = (<em>n</em>+1)/2;
        /// then for: <em>k</em> = 1, ..., <em>m</em>-1 <br/>
        /// the real part of <em>k</em>-th complex FFT coefficients is <em>x</em>[2*<em>k</em>-1];<br/>
        /// the imaginary part of <em>k</em>-th complex FFT coefficients is <em>x</em>[2*<em>k</em>].</para>
        /// <para>If <em>n</em> is even, the real of part of (<em>n</em>/2)-th complex FFT coefficients is <em>x</em>[<em>n</em>]; its imaginary part is 0.</para>
        /// The remaining complex FFT coefficients can be obtained by the symmetry relation:
        /// the (<em>n</em>-<em>k</em>)-th complex FFT coefficient is the conjugate of <em>n</em>-th complex FFT coefficient.
        /// </param>
        public void ft(double[] x)
        {
            if (x.Length != Ndim)
                throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");

            rfftf(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Forward real FFT transform. It computes the discrete transform of a real data sequence.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// <para>After FFT, <em>x</em> contains the transform coefficients used to construct <em>n</em> complex FFT coefficients.</para></param>
        /// <param name="y">The first complex (<em>n</em>+1)/2 (when <em>n</em> is odd) or (<em>n</em>/2+1) (when <em>n</em> is even) FFT coefficients.
        /// <para>The remaining complex FFT coefficients can be obtained by the symmetry relation:
        /// the (<em>n</em>-<em>k</em>)-th complex FFT coefficient is the conjugate of <em>n</em>-th complex FFT coefficient.</para></param>
        public void ft(double[] x, out Complex[] y)
        {
            if (x.Length != Ndim)
                throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");

            rfftf(Ndim, x, WaveTable);

            int yDim = (Ndim % 2 == 0) ? (Ndim / 2 + 1) : ((Ndim + 1) / 2);

            y = new Complex[yDim];

            y[0].Re = x[0];
            y[0].Im = 0.0D;

            for (int i = 1; i < (Ndim + 1) / 2; i++)
            {
                y[i].Re = x[2 * i - 1];
                y[i].Im = x[2 * i];
            }

            if (Ndim % 2 == 0)
            {
                y[Ndim / 2].Re = x[Ndim - 1];
                y[Ndim / 2].Im = 0.0D;
            }
        }

        /// <summary>
        /// Backward real FFT transform. It is the unnormalized inverse transform of <see cref="ft(double[])"/>.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed.<br/>
        /// After FFT, <em>x</em> contains the transform coefficients.</param>
        /// <remarks>See the comments of <see cref="ft(double[])"/> for the relation between <em>x</em> and complex FFT coefficients.</remarks>
        public void bt(double[] x)
        {
            if (x.Length != Ndim)
                throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");

            rfftb(Ndim, x, WaveTable);
        }

        /// <summary>
        /// Backward real FFT transform. It is the unnormalized inverse transform of <see cref="ft(double[], Complex1D)"/>.
        /// </summary>
        /// <param name="x">An array which contains the sequence to be transformed. When <em>n</em> is odd, it contains the first (<em>n</em>+1)/2 complex data; when <em>n</em> is even, it contains (<em>n</em>/2+1) complex data.</param>
        /// <param name="y">The real FFT coefficients.</param>
        /// <remarks>See the comments of <see cref="ft(double[])"/> for the relation between <em>x</em> and complex FFT coefficients.</remarks>
        public void bt(Complex[] x, out double[] y)
        {
            if (Ndim % 2 == 0)
            {
                if (x.Length != (Ndim / 2 + 1))
                    throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");
            }
            else
            {
                if (x.Length != (Ndim + 1) / 2)
                    throw new ArgumentOutOfRangeException("x", "The length of data can not match that of the wavetable");
            }

            y = new double[Ndim];
            y[0] = x[0].Re;

            for (int i = 1; i < (Ndim + 1) / 2; i++)
            {
                y[2 * i - 1] = x[i].Re;
                y[2 * i] = x[i].Im;
            }

            if (Ndim % 2 == 0)
            {
                y[Ndim - 1] = x[Ndim / 2].Re;
            }

            rfftb(Ndim, y, WaveTable);
        }
    }
}
