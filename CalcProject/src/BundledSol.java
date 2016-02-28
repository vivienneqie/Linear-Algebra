/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author nicol
 */
public class BundledSol {
    
    public Matrix m;
    public int N;
    public double[][] orig;
    
    public BundledSol(Matrix m, int N, double[][] orig) {
        this.m = m;
        this.N = N;
        this.orig = orig;
    }
    
    public String toString() {
        String toreturn = "";
        toreturn += ("Approximate Solution: \n");
        toreturn += m.toString();
        toreturn += ("\n Number of Iterations: " + N);
        toreturn += ("\n");
        return toreturn;
    }
}
