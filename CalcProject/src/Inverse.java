import java.util.Scanner;

public class Inverse {

    public static double[][] invert(double[][] arr) {
        return invert(new Matrix(arr));
    }

    public static double[][] invert(Matrix a_matrix) {

        double a[][] = a_matrix.getArray();
        int n = a.length;
        int index[] = new int[n];
        Matrix b_matrix = LinearAlgebra.identityMatrix(n);
        double b[][] = b_matrix.getArray();
        double x[][] = new double[n][n];

        upper_triangle(a, index);

        // Update the matrix b[i][j] with the ratios stored
        for (int r = 0; r < n-1; r++) {
            for (int c = r + 1; c < n; c++) {
                for (int k = 0; k < n; k++) {
                    b[index[c]][k] -= a[index[c]][r]*b[index[r]][k];
                }
            }
        }

        // Perform backward substitutions
        for (int i=0; i<n; ++i)
        {
            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
            for (int j=n-2; j>=0; --j)
            {
                x[j][i] = b[index[j]][i];
                for (int k=j+1; k<n; ++k)
                {
                    x[j][i] -= a[index[j]][k]*x[k][i];
                }
                x[j][i] /= a[index[j]][j];
            }
        }
        return x;
    }



    // partial-pivots stored in index[]
    public static void upper_triangle(double a[][], int index[]) {

        int n = index.length;
        index = fill(index, 1);
        double[] c = new double[n];

        // Find the rescaling factors, one from each row
        for (int i = 0; i < n; i++) {
            double c1 = 0;
            for (int j = 0; j < n; j++) {
                double c0 = Math.abs(a[i][j]);
                if (c0 > c1) c1 = c0;
            }
            c[i] = c1;
        }

        // Look at pivot number
        int k = 0;
        for (int j=0; j<n-1; ++j) {
            double pi1 = 0;
            for (int i=j; i<n; ++i) {
                double pi0 = Math.abs(a[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pi1) {
                    pi1 = pi0;
                    k = i;
                }
            }

            // Interchange rows according to the pivoting order
            index = swap(index, j, k);

            for (int i=j+1; i<n; ++i)
            {
                double pj = a[index[i]][j]/a[index[j]][j];

                // Record pivoting ratios below the diagonal
                a[index[i]][j] = pj;

                // Modify other elements accordingly
                for (int l = j+1; l < n; ++l)
                    a[index[i]][l] -= pj*a[index[j]][l];
            }
        }
    }

    public static int[] swap(int[] arr, int pos_i, int pos_j) {
        int temp = arr[pos_i];
        arr[pos_i] = arr[pos_j];
        arr[pos_j] = temp;
        return arr;
    }

    public static void main(String argv[])
    {
        Scanner input = new Scanner(System.in);
        System.out.println("Enter the dimension of square matrix: ");
        int n = input.nextInt();
        double a[][]= new double[n][n];
        System.out.println("Enter the elements of matrix: ");
        for(int i=0; i<n; i++)
            for(int j=0; j<n; j++)
                a[i][j] = input.nextDouble();

        double d[][] = invert(a);

        System.out.println("The inverse is: ");
        for (int i=0; i<n; ++i)
        {
            for (int j=0; j<n; ++j)
            {
                System.out.print(d[i][j]+"  ");
            }
            System.out.println();
        }
        input.close();
    }

    private static int[] fill(int[] arr, int num) {
        for (int i = 0; i < arr.length; i++) {
            arr[i] = i;
        }
        return arr;
    }
}
