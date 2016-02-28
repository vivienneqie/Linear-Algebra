import java.util.Scanner;

public class solve_lu_b {

    private Matrix L;
    private Matrix U;
    private Vector sol;
    private double error;

    public solve_lu_b(Matrix a, Vector b) {
        lu_fact LU = new lu_fact(a);
        L = LU.getL();
        U = LU.getU();
        error = LU.getError();

        double[] y = new double[b.getLength()];
        for (int i = 0; i < b.getLength(); i++) {
            y[i] = b.get(i);
            for (int j = 0; j < i; j++) {
                y[i] = y[i] - L.get(i, j) * y[j];
            }
            y[i] = y[i] / L.get(i, i);
        }

        double[] x = new double[b.getLength()];
        for (int i = b.getLength() - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < b.getLength(); j++) {
                x[i] = x[i] - U.get(i, j) * x[j];
            }
            x[i] = x[i] / U.get(i, i);
        }

        sol = new Vector(x);

    }

    public Vector getSol() {
        return sol;
    }

    public Matrix getL() {
        return L;
    }

    public Matrix getU() {
        return U;
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

        System.out.println("Enter the length of the vector: ");
        int length = input.nextInt();
        double[] v = new double[length];
        System.out.println("Enter the elements of vector V: ");
        for (int i = 0; i < length; i++) {
            v[i] = input.nextDouble();
        }

        Vector bVector = new Vector(v);

        solve_lu_b LU = new solve_lu_b(aMatrix, bVector);

        System.out.println("L:");
        System.out.println(LU.getL());
        System.out.println("U:");
        System.out.println(LU.getU());
        System.out.println("x_sol:");
        System.out.println(LU.getSol());
    }
}
