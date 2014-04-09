package barna.flux.capacitor.lp;

import barna.flux.capacitor.reconstruction.FluxCapacitor;
import lpsolve.LpSolve;
import lpsolve.LpSolveException;

import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 11/2/13
 * Time: 9:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class LPexplorer {

    public static void main(String[] args) {
        long t0= System.currentTimeMillis();
        FluxCapacitor.loadLibraries();
        LpSolve lp= createModel(2000000, 100);
        long t1= System.currentTimeMillis();
        System.out.println("time "+ (t1- t0)/ 1000+ " sec.");
        System.currentTimeMillis();
    }

    private static LpSolve createModel(int nrConstraints, int nrRestrictions) {

        try {
            LpSolve lpSolve= LpSolve.makeLp(0, nrConstraints);
            int[] idx= new int[nrConstraints];
            Random r= new Random();
            for (int i = 0; i < nrConstraints; i++) {
                idx[i]= i+ 1;
            }
            for (int i = 0; i < nrRestrictions; i++) {

                double[] val= new double[nrConstraints];
                for (int j = 0; j < val.length; j++) {
                    val[j] = r.nextDouble();
                }
                lpSolve.addConstraintex(idx.length, val, idx, LpSolve.EQ, r.nextDouble());
            }

            return lpSolve;

        } catch (LpSolveException e) {
            e.printStackTrace();
        }

        return null;
    }
}
