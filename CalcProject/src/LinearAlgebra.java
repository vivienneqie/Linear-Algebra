/**
 * Class that contains useful Linear Algebra functions
 *
 * @author Nicole Hofmann
 * @version 1.0
 */

public class LinearAlgebra {

    public static Vector matrixVectorMultiply(Matrix m, Vector v){
        double[] d = new double[v.getLength()];
        int i = 0;
        for (int r = 0; r < m.getHeight(); r++) {
            for (int c = 0; c < m.getWidth(); c++) {
                d[i] += (m.get(r, c) * v.get(c));
            }
            i++;
        }
        return new Vector(d);
    }
    
    public static Matrix identityMatrix(int dim) {
        double[][] ident = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if (i == j) {
                    ident[i][j] = 1;
                } else {
                    ident[i][j] = 0;
                }
            }
        }
        return new Matrix(ident);
    }

    public static Matrix matrixAdd(Matrix m1, Matrix m2){
        double[][] d = new double[m1.getHeight()][m1.getWidth()];
        for (int r = 0; r < m1.getHeight(); r++) {
            for (int c = 0; c < m1.getWidth(); c++) {
                d[r][c] = m1.get(r, c) + m2.get(r, c);
            }
        }
        return new Matrix(d);
    }
    
    public static Matrix matrixSubtract(Matrix m1, Matrix m2){
        double[][] d = new double[m1.getHeight()][m1.getWidth()];
        for (int r = 0; r < m1.getHeight(); r++) {
            for (int c = 0; c < m1.getWidth(); c++) {
                d[r][c] = m1.get(r, c) - m2.get(r, c);
            }
        }
        return new Matrix(d);
    }

    public static double dotProduct(Vector v1, Vector v2){
        double d = 0;
        for (int i = 0; i < v1.getLength(); i++) {
            d += (v1.get(i) * v2.get(i));
        }
        return d;
    }

    public static Vector vectorAdd(Vector v1, Vector v2){
        double[] d = new double[v1.getLength()];
        for (int i = 0; i < v1.getLength(); i++) {
            d[i] = v1.get(i) + v2.get(i);
        }
        return new Vector(d);
    }
    
    public static Vector vectorSubtract(Vector v1, Vector v2){
        double[] d = new double[v1.getLength()];
        for (int i = 0; i < v1.getLength(); i++) {
            d[i] = v1.get(i) - v2.get(i);
        }
        return new Vector(d);
    }
    
    public static Vector vectorScalarMultiple(Vector v, double s) {
        double[] d = new double[v.getLength()];
        for (int i = 0; i < v.getLength(); i++) {
            d[i] = v.get(i) * s;
        }
        return new Vector(d);
    }
    
    public static Matrix matrixScalarMultiple(Matrix m1, double s){
        double[][] d = new double[m1.getHeight()][m1.getWidth()];
        for (int r = 0; r < m1.getHeight(); r++) {
            for (int c = 0; c < m1.getWidth(); c++) {
                d[r][c] = m1.get(r, c) * s;
            }
        }
        return new Matrix(d);
    }
    
    public static Matrix matrixMultiply (Matrix m1, Matrix m2) {
        double[][]d = new double[m1.getWidth()][m2.getHeight()];
        double currSum = 0;
        for (int i = 0; i < m1.getHeight(); i++) {
            for (int j = 0; j < m2.getWidth(); j++) {
                for (int k = 0; k < m2.getHeight(); k++) {
                    currSum = currSum + m1.get(i, k) * m2.get(k, j);
                }
                d[i][j] = currSum;
                currSum = 0;
            }
        }
        return new Matrix(d);
    }
    
    public static double determinant(Matrix m) {
        double det = 0;
        double[][] matrix = m.getArray();
        int order = m.getHeight();
        if (order == 1) {
            det = matrix[0][0];
        }
        else if (order == 2) {
            det = matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1];
        } else {
            det = 0;
            for(int j1 = 0; j1 < order; j1++) {
                double[][] m1 = new double[order - 1][];
                for(int k = 0; k < (order - 1); k++) {
                    m1[k] = new double[order - 1];
                } for(int i = 1; i < order; i++) {
                    int j2 = 0;
                    for(int j = 0; j < order; j++) {
                        if(j == j1)
                            continue;
                        m1[i - 1][j2] = matrix[i][j];
                        j2++;
                    }
                }
                det += Math.pow(-1.0, 1.0 + j1 + 1.0)* matrix[0][j1] * determinant(new Matrix(m1));
            }
        }
        return det;
    }
    
    public static Matrix matrixCopy(Matrix oldMatrix) {
        int width = oldMatrix.getWidth();
        int height = oldMatrix.getHeight();
        Matrix copy = new Matrix(height, width);
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < width; j++) {
                copy.setValue(i, j, oldMatrix.get(i, j));
            }
        }
        return copy;
    }

    public static Matrix matrixPower(Matrix m, int p) {
        Matrix result = m;
        for (int i = 1; i < p; i++) {
            result = matrixMultiply(result, m);
        }
        return result;
    }
}