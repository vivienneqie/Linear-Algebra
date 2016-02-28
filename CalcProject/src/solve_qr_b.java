import java.util.Scanner;

/**
 * @author Vivienne Qie
 */
public class solve_qr_b {

    private Matrix Q_givens, Q_househ;
    private Matrix R_givens, R_househ;
    private Vector sol_givens, sol_househ;
    private double error_givens, error_househ;

    public solve_qr_b(Matrix a, Vector b) {
        qr_fact_givens givensFact = new qr_fact_givens(a);
        Q_givens = givensFact.getQ();
        R_givens = givensFact.getR();
        error_givens = givensFact.getError();

        Vector d = LinearAlgebra.matrixVectorMultiply(Q_givens.transpose(), b);

        double[] x = new double[b.getLength()];
        for (int i = b.getLength() - 1; i > -1; i--) {
            x[i] = d.get(i);
            for (int j = i + 1; j < b.getLength(); j++) {
                x[i] = x[i] - R_givens.get(i, j) * x[j];
            }
            x[i] = x[i] / R_givens.get(i, i);
        }

        sol_givens = new Vector(x);

        qr_fact_househ househFact = new qr_fact_househ(a);
        Q_househ = househFact.getQ();
        R_househ = househFact.getR();
        error_househ = househFact.getError();

        Vector d2 = LinearAlgebra.matrixVectorMultiply(Q_househ.transpose(), b);

        double[] x2 = new double[b.getLength()];
        for (int i = b.getLength() - 1; i > -1; i--) {
            x2[i] = d2.get(i);
            for (int j = i + 1; j < b.getLength(); j++) {
                x2[i] = x2[i] - R_househ.get(i, j) * x2[j];
            }
            x2[i] = x2[i] / R_househ.get(i, i);
        }

        sol_househ = new Vector(x);
    }

    public Matrix getQ_givens() {
        return Q_givens;
    }

    public Matrix getQ_househ() {
        return Q_househ;
    }

    public Matrix getR_givens() {
        return R_givens;
    }

    public Matrix getR_househ() {
        return R_househ;
    }

    public Vector getSol_givens() {
        return sol_givens;
    }

    public Vector getSol_househ() {
        return sol_househ;
    }

    public double getError_givens() {
        return error_givens;
    }

    public double getError_househ() {
        return error_househ;
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

        System.out.println("Enter the length of the vector: ");
        int length = input.nextInt();
        double[] v = new double[length];
        System.out.println("Enter the elements of vector V: ");
        for (int i = 0; i < length; i++) {
            v[i] = input.nextDouble();
        }

        Vector bVector = new Vector(v);

        solve_qr_b given = new solve_qr_b(aMatrix, bVector);

        System.out.println("Q (givens):");
        System.out.println(given.getQ_givens());
        System.out.println("R (givens):");
        System.out.println(given.getR_givens());
        System.out.println("x_sol (givens):");
        System.out.println(given.getSol_givens());
        System.out.println("Q (householder):");
        System.out.println(given.getQ_househ());
        System.out.println("R (householder):");
        System.out.println(given.getR_househ());
        System.out.println("x_sol (householder):");
        System.out.println(given.getSol_househ());
    }
}
