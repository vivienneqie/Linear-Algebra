import java.io.File;
import java.io.PrintWriter;
import java.util.*;

/**
 * @author Vivienne Qie
 */
public class power_method {

    private Vector eigenvector;
    private double eigenvalue;
    private int iteration = 0;

    public power_method(Matrix a, Vector v, double tol, int max) {

        ArrayList<Double> eigenvalues = new ArrayList<>();
        eigenvalues.add(0.0);

        while (iteration < max) {

            Matrix powerA = LinearAlgebra.matrixPower(a, iteration + 1);

            Vector vectorAu = LinearAlgebra.matrixVectorMultiply(powerA, v);
            for (int j = 0; j < vectorAu.getLength(); j++) {
                double after = vectorAu.get(j) / vectorAu.norm();
                vectorAu.setValue(j, after);
            }
            eigenvector = vectorAu;

            Vector vectorAx = LinearAlgebra.matrixVectorMultiply(a, eigenvector);
            double axx = LinearAlgebra.dotProduct(vectorAx, eigenvector);
            double xx = LinearAlgebra.dotProduct(eigenvector, eigenvector);
            eigenvalue = axx / xx;

            eigenvalues.add(eigenvalue);

            if (Math.abs(eigenvalues.get(iteration + 1)
                    - eigenvalues.get(iteration)) <= tol) {
                break;
            }

            iteration++;
        }
    }

    public Vector getEigenvector() {
        return eigenvector;
    }

    public double getEigenvalue() {
        return eigenvalue;
    }

    public int getIteration() {
        return iteration;
    }

    public static double randomWithRange(double min, double max) {
        double range = Math.abs(max - min);
        return (Math.random() * range) + (min <= max ? min : max);
    }

    public static ArrayList<Matrix> thousandMatrices() {
        ArrayList<Matrix> matrices = new ArrayList<>();
        for (int i = 0; i <  1000; i++) {
            double[][] matrix = new double[][]{{randomWithRange(-2.0, 2.0),
                    randomWithRange(-2.0, 2.0)}, {randomWithRange(-2.0, 2.0),
                    randomWithRange(-2.0, 2.0)}};
            Matrix entry = new Matrix(matrix);
            matrices.add(entry);
        }
        return matrices;
    }

    public static ArrayList<Matrix> calculateInverse(ArrayList<Matrix> matrices) {
        ArrayList<Matrix> inverses = new ArrayList<>();
        for (Matrix matrix : matrices) {
            Matrix inverse = matrix.inverse2x2();
            inverses.add(inverse);
        }
        return inverses;
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

        Vector vVector = new Vector(v);

        System.out.println("Enter error tolerance for power method: ");
        double tolerance = input.nextDouble();

        System.out.println("Enter max number of iterations before quitting: ");
        int max_iter = input.nextInt();

        power_method pm = new power_method(aMatrix, vVector, tolerance, max_iter);

        System.out.println("Dominant eigenvalue: ");
        System.out.println(pm.getEigenvalue());
        System.out.println("Dominant eigenvector: ");
        System.out.println(pm.getEigenvector());

        if (pm.getIteration() >= max_iter) {
            System.out.println("Did not converge after " + max_iter + " iterations.");
        } else {
            System.out.println("Took " + pm.getIteration() + " iterations to converge.");
        }
    }
}