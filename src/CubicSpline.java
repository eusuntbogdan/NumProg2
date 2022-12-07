import java.util.Arrays;

/**
 * Die Klasse CubicSpline bietet eine Implementierung der kubischen Splines. Sie
 * dient uns zur effizienten Interpolation von aequidistanten Stuetzpunkten.
 *
 * @author braeckle
 */
public class CubicSpline implements InterpolationMethod {

    /**
     * linke und rechte Intervallgrenze x[0] bzw. x[n]
     */
    double a, b;

    /**
     * Anzahl an Intervallen
     */
    int n;

    /**
     * Intervallbreite
     */
    double h;

    /**
     * Stuetzwerte an den aequidistanten Stuetzstellen
     */
    double[] y;

    /**
     * zu berechnende Ableitunge an den Stuetzstellen
     */
    double yprime[];

    /**
     * {@inheritDoc} Zusaetzlich werden die Ableitungen der stueckweisen
     * Polynome an den Stuetzstellen berechnet. Als Randbedingungen setzten wir
     * die Ableitungen an den Stellen x[0] und x[n] = 0.
     */
    @Override
    public void init(double a, double b, int n, double[] y) {
        this.a = a;
        this.b = b;
        this.n = n;
        h = ((double) b - a) / (n);

        this.y = Arrays.copyOf(y, n + 1);

        /* Randbedingungen setzten */
        yprime = new double[n + 1];
        yprime[0] = 0;
        yprime[n] = 0;

        /* Ableitungen berechnen. Nur noetig, wenn n > 1 */
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * getDerivatives gibt die Ableitungen yprime zurueck
     */
    public double[] getDerivatives() {
        return yprime;
    }

    /**
     * Setzt die Ableitungen an den Raendern x[0] und x[n] neu auf yprime0 bzw.
     * yprimen. Anschliessend werden alle Ableitungen aktualisiert.
     */
    public void setBoundaryConditions(double yprime0, double yprimen) {
        yprime[0] = yprime0;
        yprime[n] = yprimen;
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * Berechnet die Ableitungen der stueckweisen kubischen Polynome an den
     * einzelnen Stuetzstellen. Dazu wird ein lineares System Ax=c mit einer
     * Tridiagonalen Matrix A und der rechten Seite c aufgebaut und geloest.
     * Anschliessend sind die berechneten Ableitungen y1' bis yn-1' in der
     * Membervariable yprime gespeichert.
     * <p>
     * Zum Zeitpunkt des Aufrufs stehen die Randbedingungen in yprime[0] und yprime[n].
     * Speziell bei den "kleinen" Faellen mit Intervallzahlen n = 2
     * oder 3 muss auf die Struktur des Gleichungssystems geachtet werden. Der
     * Fall n = 1 wird hier nicht beachtet, da dann keine weiteren Ableitungen
     * berechnet werden muessen.
     */
    public void computeDerivatives() {
        double[] c = new double[n - 1];
        double[] a = new double[n - 1];
        double[] b = new double[n - 1];
        double[] d = new double[n - 1];

        //Compute solution vector
        for (int i = 1; i < n - 2; i++) {
            d[i] = y[i + 2] - y[i];
        }
        for (int i = 1; i < n - 2; i++) {
            d[i] *= 3 / h;
        }
        d[0] = (3/h) * (y[2] - y[0] - (h/3) * yprime[0]);
        d[n-2] = (3/h) * (y[n] - y[n-2] - (h/3) * yprime[n]);

        //Initialise tridiagonal matrix
        for (int i = 0; i < n - 1; i++) {
            c[i] = 1;
            b[i] = 4;
            a[i] = 1;
        }

        //Solve tridiagonal matrix equation A * v = x
        //Result will be in x at the end
        for (int i = 1; i < d.length; i++) {
            double m = a[i - 1] / b[i-1];
            b[i] = b[i] - m * c[i-1];
            d[i] = d[i] - d[i-1] * m;
        }

        double[] x = Arrays.copyOf(b, b.length);
        x[x.length - 1] = d[d.length - 1] / b[b.length - 1];

        for (int i = d.length - 2; i >= 0; i--) {
            x[i] = (d[i] - c[i] * x[i+1])/b[i];
        }
        System.arraycopy(x, 0, yprime, 1, x.length);
    }


    /**
     * {@inheritDoc} Liegt z ausserhalb der Stuetzgrenzen, werden die
     * aeussersten Werte y[0] bzw. y[n] zurueckgegeben. Liegt z zwischen den
     * Stuetzstellen x_i und x_i+1, wird z in das Intervall [0,1] transformiert
     * und das entsprechende kubische Hermite-Polynom ausgewertet.
     */
    @Override
    public double evaluate(double z) {

        if (z <= a) {
            return y[0];
        }
        if (z >= b) {
            return y[n];
        }
        int i;
        for (i = 0; i < n; i++) {
            if (z >= a + i * h && z <= a + (i + 1) * h) {
                break;
            }
        }
        double t = (z - (a + i * h)) / h;
        double result0 = y[i] * (1 + (-3) * t * t + 2 * t * t * t);
        double result1 = y[i + 1] * (3 * t * t - 2 * t * t * t);
        double result2 = h * yprime[i] * (t - 2 * t * t + t * t * t);
        double result3 = h * yprime[i + 1] * (-t * t + t * t * t);
        return result0 + result1 + result2 + result3;
    }
}
