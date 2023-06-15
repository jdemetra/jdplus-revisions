package jdplus.revisions.base.core.treatment;

import java.util.ArrayList;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.api.stats.StatException;

/**
 *
 * @author clemasso
 */

public class PreTreatment {
   
    private DoubleSeq a; 
    private DoubleSeq b;
    private Matrix m;
    
    public PreTreatment (DoubleSeq a, DoubleSeq b){
        this.a = a;
        this.b = b;
    }
    
    public PreTreatment (Matrix m){
        this.m = m;
    }
    
    public DoubleSeq getDoubleSeqa(){
        return this.a;
    }
    
    public DoubleSeq getDoubleSeqb(){
        return this.b;
    }
    
    public Matrix getMatrix (){
        return this.m;
    }
    
    /**
     * The method compares two sequences by element and  
     * removes the corresponding element in each sequence if at 
     * least one of them is missing
     *
     * @param a The first sequence
     * @param b The second sequence
     * @return matrix including the two cleaned sequences 
     */
    
    public static Matrix cleanNaN(DoubleSeq a, DoubleSeq b){
        
        int n = a.length();
        if (b.length() != n) {
            throw new StatException("Non compatible data");
        }
        
        double[] aAr = a.toArray();
        double[] bAr = b.toArray();

        ArrayList<Double> aAL= new ArrayList<>();
        ArrayList<Double> bAL = new ArrayList<>();

        for(int k = 0; k < n; ++k){
            if(Double.isFinite(aAr[k]) && Double.isFinite(bAr[k])){  
                aAL.add(aAr[k]);
                bAL.add(bAr[k]);
            }   
        }
        
        ArrayList<Double> cAL = new ArrayList<>();
        cAL.addAll(aAL);
        cAL.addAll(bAL);
        double[] cAr = cAL.stream().mapToDouble(d -> d).toArray(); 
        
        return(Matrix.of(cAr, aAL.size(), 2));
    }
    
    /**
     * The method removes missing value from a sequence
     *
     * @param a a sequence
     * @return a DoubleSeq object with the cleaned series 
     */
    
    public static DoubleSeq cleanNaN(DoubleSeq a){
        
        int n = a.length();
        double[] aAr = a.toArray();
        ArrayList<Double> aAL= new ArrayList<>();

        for(int k = 0; k < n; ++k){
            if(Double.isFinite(aAr[k])){  
                aAL.add(aAr[k]);
            }   
        }
        double[] cAr = aAL.stream().mapToDouble(d -> d).toArray(); 
        
        return(DoubleSeq.of(cAr));
    }
    
    
       /**
     * The method compares x sequences by element and  
     * removes the corresponding element in each sequence if at 
     * least one of them is missing
     *
     * @param m Matrix of sequences
     * @return matrix including the two cleaned sequences 
     */
    
    public static Matrix cleanNaN(Matrix m){
        
        int nr = m.getRowsCount();
        int nc = m.getColumnsCount();

        double[][] mAr = new double[nr][nc];

        for (int i = 0; i < nr; ++i) {
            for (int j = 0; j < nc; ++j) {
                mAr[i][j] = m.get(i, j);
            }
        }

        ArrayList<Double> mAL = new ArrayList<>();
        boolean keepRow;
        int nrc = 0;
        for (int i = 0; i < nr; ++i) {
            keepRow = true;
            for (int j = 0; j < nc; ++j) {
                if (Double.isNaN(mAr[i][j])) {
                    keepRow = false;
                    break;
                }
            }
            if (keepRow) {
                for (int j = 0; j < nc; ++j) {
                    mAL.add(mAr[i][j]);
                }
                ++nrc;
            }
        }

        double[] mcAr = mAL.stream().mapToDouble(d -> d).toArray();
        
        // Re-arranging array
        int l=mcAr.length;
        double[] mcoAr = new double[l];
        for (int i = 0; i<l; ++i) {
            int b = i/nc;
            mcoAr[(i*nrc)-(b*(l-1))]=mcAr[i];
        } 

        return(Matrix.of(mcoAr, nrc, nc));
    }
    
    /**
     * The method replaces all missing values in a matrix by  
     * a the byValue 
     *
     * @param m a matrix
     * @param byValue value used to replace the missing ones 
     * @return cleaned matrix
     */
    
    public static Matrix cleanNaN(Matrix m, double byValue){
        
        double[] mAr = m.toArray();
        int n = mAr.length;
        
        for (int k = 0; k < n; ++k) {
            if (Double.isNaN(mAr[k])) {
                mAr[k] = byValue;
            }
        }    
        return(Matrix.of(mAr, m.getRowsCount(), m.getColumnsCount()));
    }
}