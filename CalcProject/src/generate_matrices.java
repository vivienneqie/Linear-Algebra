import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 * @author Vivienne Qie
 */
public class generate_matrices {

    public static void main(String[] args) throws Exception {
        ArrayList<Matrix> matrices = power_method.thousandMatrices();
        ArrayList<Matrix> inverses = power_method.calculateInverse(matrices);

        double[] vector = new double[]{1, 0};
        Vector v = new Vector(vector);
        double tol = 0.00005;
        int n = 100;

        power_method regPowerMethod;
        power_method inversePowerMethod;

        ArrayList<Double> largestEigenvalue = new ArrayList<>();
        ArrayList<Double> smallestEigenvalue = new ArrayList<>();
        ArrayList<Integer> iterA = new ArrayList<>();
        ArrayList<Integer> iterInverse = new ArrayList<>();
        ArrayList<Double> determinant = new ArrayList<>();
        ArrayList<Double> trace = new ArrayList<>();

        for (Matrix m : matrices) {
            regPowerMethod = new power_method(m, v, tol, n);
            largestEigenvalue.add(regPowerMethod.getEigenvalue());
            iterA.add(regPowerMethod.getIteration());
            determinant.add(LinearAlgebra.determinant(m));
            trace.add(m.trace());
        }

        for (Matrix m : inverses) {
            inversePowerMethod = new power_method(m, v, tol, n);
            smallestEigenvalue.add(1 / inversePowerMethod.getEigenvalue());
            iterInverse.add(inversePowerMethod.getIteration());
        }

        PrintWriter largeEig = new PrintWriter(new File("largestEig.txt"));
        for (double i : largestEigenvalue) {
            largeEig.println(i);
        }
        largeEig.close();

        PrintWriter smallEig = new PrintWriter(new File("smallestEig.txt"));
        for (double i : smallestEigenvalue) {
            smallEig.println(i);
        }
        smallEig.close();

        PrintWriter determinantFile = new PrintWriter(new File("determinants.txt"));
        PrintWriter traceFile = new PrintWriter(new File("traces.txt"));
        PrintWriter fileAIter = new PrintWriter(new File("iterations_of_A.txt"));
        PrintWriter fileInverse = new PrintWriter(new File("iterations_of_A_inverse.txt"));

        for (double i : determinant) {
            determinantFile.println(i);
        }
        for (double j : trace) {
            traceFile.println(j);
        }
        for (int k : iterA) {
            fileAIter.println(k);
        }
        for (int l : iterInverse) {
            fileInverse.println(l);
        }

        determinantFile.close();
        traceFile.close();
        fileAIter.close();
        fileInverse.close();

        PrintWriter matrix = new PrintWriter(new File("A.txt"));
        for (Matrix m : matrices) {
            matrix.println(m.toString());
            matrix.println();
        }
        matrix.close();

        PrintWriter inverse = new PrintWriter(new File("inverses.txt"));
        for (Matrix m : inverses) {
            inverse.println(m.toString());
            inverse.println();
        }
        inverse.close();

    }
}
