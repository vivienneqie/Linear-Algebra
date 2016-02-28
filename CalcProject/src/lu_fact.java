import java.util.Scanner;

public class lu_fact {

    private Matrix L, U;
    private double error;

    public lu_fact(Matrix a) {
        int rows = a.getHeight();
        int columns = a.getWidth();

        L = LinearAlgebra.identityMatrix(rows);
        U = a;

        for (int i = 0; i < rows - 1; i++) {
            for (int j = i + 1; j < rows; j++) {
                L.setValue(j, i, U.get(j, i) / U.get(i, i));
                for (int k = i; k < rows; k++) {
                    U.setValue(j, k, U.get(j, k) - (L.get(j, i) * U.get(i, k)));
                }
            }
        }
        Matrix LU = LinearAlgebra.matrixMultiply(L, U);
        Matrix LU_A = LinearAlgebra.matrixSubtract(LU, a);

        double max = 0.0;
        for (int i = 0; i < LU_A.getHeight(); i++) {
            for (int j = 0; j < LU_A.getWidth(); j++) {
                if (Math.abs(LU_A.get(i, j)) > max) {
                    max = Math.abs(LU_A.get(i, j));
                }
            }
        }
        error = max;
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

        lu_fact lu = new lu_fact(aMatrix);

        System.out.println("L:");
        System.out.println(lu.getL());
        System.out.println("U:");
        System.out.println(lu.getU());
        System.out.println("Error: " + lu.getError());
    }

}