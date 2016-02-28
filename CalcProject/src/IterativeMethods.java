import java.util.ArrayList;
import java.util.Scanner;

/**
 *
 * @author Nicole
 */
public class IterativeMethods {
    //private static final double[][] A = new double[][]{ {1, (1/2), (1/3)},
            //{(1/2), (1), (1/4)}, {(1/3), (1/4), 1} };
    private static double[][] A;
    //private static final double[] b = new double[]{ 0.1, 0.1, 0.1 };
    private static double[] b;
    private static double tol = 0.00005;
    private static int M1 = 100;
    private static final double[] xexact = new double[]{ (9/190), (28/475), (33/475) };
    private static int NJacobiCurr;
    private static int NGSCurr;
    private static double NJacobiTotal = 0;
    private static double NGSTotal = 0;
    
    public static void main(String[] args) {
        //
        Scanner input = new Scanner(System.in);
        System.out.println("Enter the dimension of the square matrix: ");
        int n = input.nextInt();
        A = new double[n][n];
        System.out.println("Enter the elements of matrix A: ");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = input.nextDouble();
            }
        }

        System.out.println("Enter the length of the vector b: ");
        int length = input.nextInt();
        b = new double[length];
        System.out.println("Enter the elements of vector b: ");
        for (int i = 0; i < length; i++) {
            b[i] = input.nextDouble();
        }

        System.out.println("Enter error tolerance for power method: ");
        tol = input.nextDouble();

        System.out.println("Enter max number of iterations before quitting: ");
        M1 = input.nextInt();
        
        //
        double[] jacobiX;
        double[] gsX;
        ArrayList<double[]> matrixList = generateMatricies(100);
        ArrayList<double[]> jacobiResList = new ArrayList<>();
        ArrayList<double[]> gsResList = new ArrayList<>();
        
        for (int i = 0; i < 100; i++) {
            jacobiX = jacobi_iter(matrixList.get(i), tol, M1);
            jacobiResList.add(jacobiX);
            NJacobiTotal += NJacobiCurr;
            gsX = gs_iter(matrixList.get(i), tol, M1);
            gsResList.add(gsX);
            NGSTotal += NGSCurr;
        }
        
        double[] jacobiApproxSol = approxSol(jacobiResList);
        double[] gsApproxSol = approxSol(gsResList);
        
        System.out.println("");
        System.out.println("Jacobi Approximation: ");
        System.out.println(new Vector(jacobiApproxSol).toString());
        System.out.println("Gauss-Seidel Approximation");
        System.out.println(new Vector(gsApproxSol).toString());
        System.out.println("");
        
        double jacobiError = calcError(jacobiApproxSol);
        double gsError = calcError(gsApproxSol);
        
        System.out.println("Jacobi Error: " + jacobiError);
        System.out.println("G-S Error: " + gsError);
        
        NJacobiTotal /= 100;
        NGSTotal /= 100;
        double NRatio = NJacobiTotal/NGSTotal;
        
        System.out.println("Ratio of Jacobi Iterations to G-S Iterations: " + NRatio);
    }
    
    private static double calcError(double[] approx) {
        double[] temp = vectorSubtract(approx, xexact);
        Vector t = new Vector(temp);
        return t.norm();
    }
    
    private static double[] approxSol(ArrayList<double[]> resList) {
        double[] result = new double[3];
        for (int i = 0; i < resList.size(); i++) {
            result = vectorAdd(result, resList.get(i));
        }
        result[0] /= 100;
        result[1] /= 100;
        result[2] /= 100;
        return result;
    }
    
    private static ArrayList<double[]> generateMatricies(int x) {
        ArrayList<double[]> mList = new ArrayList<>();
        for (int i = 0; i < x; i++) {
            double[] d = new double[3];
            d[0] = (Math.random()*2 - 1);
            d[1] = (Math.random()*2 - 1);
            d[2] = (Math.random()*2 - 1);
            mList.add(d);
        }
        return mList;
    }
    
    private static double[] jacobi_iter(double[] x0, double tolerance, int M) {
        double[] x;
        double[][] L = new double[A.length][A.length];
        double[][] U = new double[A.length][A.length];
        double[][] D = new double[A.length][A.length];
        makeLUD(L, U, D);
        int iters = 0;
        
        double[][] D_inv = jacobi_inverse(D);
        double[][] LU = jacobi_LU(L, U);
        
        double[] k = x0;
        double[] k1 = new double[3];
        
        while (iters < M) {
            double[] temp = mvMult(LU, k);
            temp = vectorAdd(temp, b);
            temp = mvMult(D_inv, temp);
            k1 = temp;
            iters++;
            
            //check w/ tolerance if you should stop
            double[] temp2 = vectorSubtract(k1, k);
            Vector t = new Vector(temp2);
            double error = t.norm();
            
            
            if (error <= tolerance) {
                break;
            }
            k = k1;
            
            NJacobiCurr = iters;
        }

        x = k1;
        return x;
    }
    
    private static double[] gs_iter(double[] x0, double tolerance, int M) {
        double[] x = new double[3];
        double[][] L = new double[A.length][A.length];
        double[][] U = new double[A.length][A.length];
        double[][] D = new double[A.length][A.length];
        makeLUD(L, U, D);
        int iters = 0;
        
        double[][] LD_inv = gs_inverse(L, D);
        double[][] U_neg = gs_U(U);
        
        double[] k = x0;
        double[] k1 = new double[3];
        
        while (iters < M) {
            double[] temp = mvMult(U_neg, k);
            temp = vectorAdd(temp, b);
            temp = mvMult(LD_inv, temp);
            k1 = temp;
            iters++;
            
            //check w/ tolerance if you should stop
            double[] temp2 = vectorSubtract(k1, k);
            Vector t = new Vector(temp2);
            double error = t.norm();
           
            if (error <= tolerance) {
                break;
            }
            k = k1;
        }
        
        NGSCurr = iters;

        x = k1;
        return x;
    }
    
    private static double[][] jacobi_inverse(double[][] D) {
        double[][] result = new double[3][3];
        for (int r = 0; r < D.length; r++) {
            for (int c = 0; c < D[0].length; c++) {
                if (r == c) {
                    if (D[r][c] == 0) {
                        result[r][c] = 0;
                    } else {
                        result[r][c] = 1/(D[r][c]);
                    }
                } else {
                    result[r][c] = 0;
                }
            }
        }
        return result;
    }
    
    private static double[][] jacobi_LU(double[][] L, double[][] U) {
        double[][] result = new double[3][3];
        // -L-U
        for (int r = 0; r < result.length; r++) {
            for (int c = 0; c < result.length; c++) {
                result[r][c] = ((-1)*(L[r][c])) - U[r][c];
            }
        }
        return result;
    }
    
    private static double[][] gs_inverse(double[][] L, double[][] D) {
        double[][] temp = new double[3][3];
        double[][] result = new double[3][3];
        //Add L + D
        for (int r = 0; r < L.length; r++) {
            for (int c = 0; c < L.length; c++) {
                temp[r][c] = L[r][c] + D[r][c];
            }
        }
        //Use Formula
        result[0][0] = 1/temp[0][0];
        result[0][1] = 0;
        result[0][2] = 0;
        result[1][0] = (-1*temp[1][0])/(temp[0][0]*temp[1][1]);
        result[1][1] = 1/temp[1][1];
        result[1][2] = 0;
        result[2][0] = ((-1*temp[1][1]*temp[2][0]) + (temp[1][0]*temp[2][1]))/(temp[0][0]*temp[1][1]*temp[2][2]);
        result[2][1] = (-1*temp[2][1])/(temp[1][1]*temp[2][2]);
        result[2][2] = 1/temp[2][2];
        return result;
    }
    
    private static double[][] gs_U(double[][] U) {
        double[][] result = new double[3][3];
        for (int r = 0; r < U.length; r++) {
            for (int c = 0; c < U.length; c++) {
                result[r][c] = (-1*U[r][c]);
            }
        }
        return result;
    }
    
    private static void makeLUD(double[][] L, double[][] U, double[][] D) {
        //Set up D
        for (int r = 0; r < A.length; r++) {
            for (int c = 0; c < A[0].length; c++) {
                if (r == c) {
                    D[r][c] = A[r][c];
                } else {
                    D[r][c] = 0;
                }
            }
        }
        //Set up L
        L[0][0] = 0;
        L[0][1] = 0;
        L[0][2] = 0;
        L[1][0] = A[1][0];
        L[1][1] = 0;
        L[1][2] = 0;
        L[2][0] = A[2][0];
        L[2][1] = A[2][1];
        L[2][2] = 0;
        //Set up U
        U[0][0] = 0;
        U[0][1] = A[0][1];
        U[0][2] = A[0][2];
        U[1][0] = 0;
        U[1][1] = 0;
        U[1][2] = A[1][2];
        U[2][0] = 0;
        U[2][1] = 0;
        U[2][2] = 0;
    } 
    
    public static double[] mvMult(double[][] m, double[] v) {
        double[] result = new double[3];
        int i = 0;
        for (int r = 0; r < m.length; r++) {
            for (int c = 0; c < m[0].length; c++) {
                result[i] += (m[r][c] * v[c]);
            }
            i++;
        }
        return result;
    }
    
    public static double[] vectorAdd(double[] v1, double[] v2) {
        double[] result = new double[3];
        for (int i = 0; i < v1.length; i++) {
            result[i] = v1[i] + v2[i];
        }
        return result;
    }

    public static double[] vectorSubtract(double[] v1, double[] v2) {
        double[] result = new double[3];
        for (int i = 0; i < v1.length; i++) {
            result[i] = v1[i] - v2[i];
        }
        return result;
    }
}