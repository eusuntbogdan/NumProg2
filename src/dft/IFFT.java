package dft;

import java.util.Arrays;

/**
 * Schnelle inverse Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class IFFT {
    /**
     * Schnelle inverse Fourier-Transformation (IFFT).
     *
     * Die Funktion nimmt an, dass die Laenge des Arrays c immer eine
     * Zweierpotenz ist. Es gilt also: c.length == 2^m fuer ein beliebiges m.
     */
    public static Complex[] ifft(Complex[] c) {
        Complex[] invers = new Complex[c.length];
        if(c.length == 1){
            invers[0] = c[0];
        }
        else{
            int m = c.length/2;
            // make two arrays one with all numbers ind even places of c and one with all numbers in uneven places of c
            Complex[] evenPart = new Complex[m];
            Complex[] unevenPart = new Complex[m];
            for(int i = 0; i < m; i++){
                evenPart[i] = c[i*2];
                unevenPart[i] = c[i*2+1];
            }
            evenPart = ifft(evenPart);
            unevenPart = ifft(unevenPart);
            Complex omega = new Complex(Math.cos(2*Math.PI*1/c.length),Math.sin(2*Math.PI*1/c.length));
            for(int j = 0; j < m; j++ ){
                Complex secondPart = unevenPart[j].mul(omega);
                invers[j] = evenPart[j].add(secondPart);
                invers[m+j] = evenPart[j].sub(secondPart);
            }
        }
        return invers;
    }
}
