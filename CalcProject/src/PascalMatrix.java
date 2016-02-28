//import java.io.File;
//import java.io.PrintWriter;
//import java.util.*;

public class PascalMatrix {
    private int n;
    private Matrix pascal;
    public PascalMatrix(int n) {
        this.n = n;
        pascal = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            pascal.setValue(0, i, 1);
            pascal.setValue(i, 0, 1);
        }
        for (int row = 1; row < n; row++) {
            for (int col = 1; col < n; col++) {
                double newVal = (pascal.get(row-1, col))+(pascal.get(row, col-1));
                pascal.setValue(row, col, newVal);
            }
        }
    }
    public Matrix getPascal() {
        return pascal;
    }
    
    private static Vector getB(int i) {
        Vector b = new Vector(i);
        for (int k = 0; k < i; k++) {
            b.setValue(k, 1);
        }
        return b;
    }
    
    public static void main(String[] args) {
        System.out.println("We chose to show the error for LU Factorization!");
        System.out.println("Here are all the Pascal Matrices from n=2 to n=12");
        
        for (int i = 2; i <= 12; i++) {
            System.out.println("When n = " + i + ", Pascal Matrix is : ");
            Matrix input = (new PascalMatrix(i).getPascal());
            Vector b = getB(i);
            solve_lu_b solve = new solve_lu_b(input, b);
            Vector x = (solve.getSol());
//            System.out.println("Here is the solution x :");
//            System.out.println(x.toString());
//            System.out.println("\n Here is the error :");
            double err = (solve.getError());
            System.out.println(i + " " + err + "\n");
            
            System.out.println("/////////////////////");
            Vector px = LinearAlgebra.matrixVectorMultiply(input, x);
            Vector pxb = LinearAlgebra.vectorSubtract(px, b);
            System.out.println(i + " " + pxb.norm() + "\n");
            System.out.println("/////////////////////");
        }
        
    }
}