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
    /// Complex number.
    /// </summary>
    public struct Complex
    {
        /// <summary>
        /// Real part.
        /// </summary>
        public double Re;
        /// <summary>
        /// Imaginary part.
        /// </summary>
        public double Im;

        public Complex(double re, double im)
        {
            this.Re = re;
            this.Im = im;
        }

        public void Conjugate()
        {
            this.Im = -this.Im;
        }

        public Complex Conjugate(Complex complex)
        {
            return new Complex(complex.Re, -complex.Im);
        }

        public override string ToString()
        {
            return String.Format("{0}+{1}j", this.Re, this.Im);
        }
    }
}
