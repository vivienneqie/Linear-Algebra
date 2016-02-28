import java.util.ArrayList;
import java.util.Scanner;

/**
 * @author Vivienne Qie
 */
public class qr_fact_givens {

    private Matrix Q;
    private Matrix R;
    private double error;

    public qr_fact_givens(Matrix a) {
        R = a;
        ArrayList<Matrix> givensTrans = new ArrayList<>();

        for (int i = 0; i < a.getWidth(); i++) {
            for (int j = a.getHeight() - 1; j > i; j--) {
                Matrix g = LinearAlgebra.identityMatrix(a.getWidth());
                double cos = generateCos(R.get(j - 1, i), R.get(j, i));
                double sin = generateSin(R.get(j - 1, i), R.get(j, i));
                g.setValue(j, j, cos);
                g.setValue(j, j - 1, sin);
                g.setValue(j - 1, j, -sin);
                g.setValue(j - 1, j - 1, cos);
                Matrix gt = g.transpose();
                givensTrans.add(gt);
                R = LinearAlgebra.matrixMultiply(g, R);
            }
        }

        Q = LinearAlgebra.identityMatrix(a.getWidth());

        for (Matrix m : givensTrans) {
            Q = LinearAlgebra.matrixMultiply(Q, m);
        }

        Matrix QR = LinearAlgebra.matrixMultiply(Q, R);
        Matrix QR_A = LinearAlgebra.matrixSubtract(QR, a);

        double max = 0.0;
        for (int i = 0; i < QR_A.getHeight(); i++) {
            for (int j = 0; j < QR_A.getWidth(); j++) {
                if (Math.abs(a.get(i, j)) > max) {
                    max = Math.abs(a.get(i, j));
                }
            }
        }
        error = max;
    }

    public double generateCos(double a, double b) {
        double r = Math.sqrt(a * a + b * b);
        double cos = a / r;
        return cos;
    }

    public double generateSin(double a, double b) {
        double r = Math.sqrt(a * a + b * b);
        double sin = -b / r;
        return sin;
    }

    public Matrix getQ() {
        return Q;
    }

    public Matrix getR() {
        return R;
    }

    public double getError() {
        return error;
    }


    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);
        System.out.println("Enter the dimension of the square matrix: ");
        int n = input.nextInt();
        double[][] a = new double[n][n];
        System.out.println("Enter the elements of matrix A: ");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = input.nextDouble();
            }
        }

        Matrix aMatrix = new Matrix(a);

        qr_fact_givens givens = new qr_fact_givens(aMatrix);

        System.out.println("Q:");
        System.out.println(givens.getQ());
        System.out.println("R:");
        System.out.println(givens.getR());
        System.out.println("Error: " + givens.getError());
    }
}
