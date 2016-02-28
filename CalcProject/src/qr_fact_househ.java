import java.util.Scanner;

public class qr_fact_househ {

    private double[][] QR;
    private double[] diagOfR;
    private int n;
    private double error;
    private Matrix matrixA;

    public qr_fact_househ (Matrix A) {

        matrixA = A;
       
        QR = A.getArrayCopy();
        if (QR.length != QR[0].length) {
            throw new IllegalArgumentException("'A' can only be a square "
                    + "matrix");
        }
        n = QR.length;   
        diagOfR = new double[n];

        for (int c = 0; c < n; c++) {
            // The norm of c-th column 
            double norm = calculateNorm(c);
            if (norm != 0.0) {
                // Householder Vector
                if (QR[c][c] < 0) {
                   norm = -norm;
                }
                for (int i = c; i < n; i++) {
                   QR[i][c] /= norm;
                }
                QR[c][c] += 1.0;
                // Transform other columns accordingly
                transform_H_Col(c);
            }
            diagOfR[c] = -norm;
        }
    }

    public Matrix getR () {
      Matrix X = new Matrix(n,n);
      double[][] R = X.getArray();
      for (int i = 0; i < n; i++) {
         for (int j = 0; j < n; j++) {
            if (i < j) {
               R[i][j] = QR[i][j];
            } else if (i == j) {
               R[i][j] = diagOfR[i];
            } else {
               R[i][j] = 0.0;
            }
         }
      }
      return X;
   }

    public Matrix getQ () {
        Matrix X = new Matrix(n,n);
        double[][] Q = X.getArray();
        for (int k = n-1; k >= 0; k--) {
            for (int i = 0; i < n; i++) {
                Q[i][k] = 0.0;
            }
            Q[k][k] = 1.0;
            for (int j = k; j < n; j++) {
                if (QR[k][k] != 0) {
                    double s = 0.0;
                    for (int i = k; i < n; i++) {
                        s += QR[i][k]*Q[i][j];
                    }
                    s = -s/QR[k][k];
                    for (int i = k; i < n; i++) {
                        Q[i][j] += s*QR[i][k];
                    }
                }
            }   
        }
       
        return X;
    }   
   
    public Matrix getH () {
        Matrix X = new Matrix(n,n);
        double[][] H = X.getArray();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i >= j) {
                    H[i][j] = QR[i][j];
                } else {
                    H[i][j] = 0.0;
                }
            }
        }
        return X;
   }
   
    private double calculateNorm(int c) {
        double norm = 0.0;
        for (int r = c; r < n; r++) {
            norm = Math.hypot(QR[r][c], norm);   // sqrt((QR[r][c])^2+(norm)^2)
        }
        return norm;
    }
    
    private void transform_H_Col(int c) {
        for (int j = c+1; j < n; j++) {
            double factor = 0.0; 
            for (int i = c; i < n; i++) {
                factor += QR[i][c] * QR[i][j];
            }
            factor = -(factor / QR[c][c]);
            for (int i = c; i < n; i++) {
                QR[i][j] += factor * QR[i][c];
            }            
        }
    }

    public Matrix getA() {
        return matrixA;
    }
    
    public Matrix qr_fact_househ (Matrix b) {
    // check dimension
        for (int i = 0; i < n; i++) {
            if (diagOfR[i] == 0) {
                throw new RuntimeException("Not enough ranks");
            }
        }
        if (b.getHeight() != n) {
            throw new IllegalArgumentException("Row dimensions have to be the "
                    + "same");
        }
       
        double[][] x = b.getArrayCopy();     
        int b_width = b.getWidth();
        for (int c = 0; c < n; c++) {
            for (int r = 0; r < b_width; r++) {
                double factor = 0.0;
                for (int i = c; i < n; i++) {
                    factor += x[i][r] * QR[i][c];
                }
                factor = -factor/QR[c][c];
                for (int i = c; i < n; i++) {
                    x[i][r] += factor * QR[i][c]; 
                }
            }
        }
        // Solve R*X = Y;
      for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < b_width; j++) {
            x[k][j] /= diagOfR[k];
         }
         for (int i = 0; i < k; i++) {
            for (int j = 0; j < b_width; j++) {
               x[i][j] -= x[k][j]*QR[i][k];
            }
         }
      }
      return (new Matrix(x).getMatrix(0,n-1,0,b_width-1));
    }

    public double getError() {
        Matrix QR = LinearAlgebra.matrixMultiply(getQ(), getR());
        Matrix QR_A = LinearAlgebra.matrixSubtract(QR, getA());

        double max = 0.0;
        for (int i = 0; i < QR_A.getHeight(); i++) {
            for (int j = 0; j < QR_A.getWidth(); j++) {
                if (Math.abs(getA().get(i, j)) > max) {
                    max = Math.abs(getA().get(i, j));
                }
            }
        }
        error = max;

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

        qr_fact_househ househ = new qr_fact_househ(aMatrix);

        System.out.println("Q:");
        System.out.println(househ.getQ());
        System.out.println("R:");
        System.out.println(househ.getR());
        System.out.println("Error: " + househ.getError());
    }

}