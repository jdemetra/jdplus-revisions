/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.revisions.base.core.treatment;

import jdplus.toolkit.base.api.math.matrices.Matrix;
import static org.assertj.core.api.Assertions.assertThat;
import static org.assertj.core.api.Assertions.byLessThan;
import org.junit.jupiter.api.Test;
import jdplus.revisions.base.core.treatment.PreTreatment;

/**
 *
 * @author LEMASSO
 */
public class PreTreatmentTest {
    
    //private static final DoubleSeq a = DoubleSeq.of(-1,1,-1,1,-1,1,-1,1,-1,1);
    //private static final DoubleSeq b = DoubleSeq.of(1,0,0,0,0,0,0,0,0,.1);


     @Test
    public void testMatrixInput() {
        
        double[] dataInput = {1.1,1.2,1.3,1.4,1.5,1.6,Double.NaN,
                              2.1,2.2,2.3,2.4,2.5,2.6,2.7,
                              3.1,3.2,Double.NaN,3.4,3.5,Double.NaN,Double.NaN};
        Matrix mInput = Matrix.of(dataInput, 7, 3);
        double[] dataOutput = {1.1,1.2,1.4,1.5,
                               2.1,2.2,2.4,2.5,
                               3.1,3.2,3.4,3.5};
        Matrix mOutput = Matrix.of(dataOutput, 4, 3);
        
        assertThat(PreTreatment.cleanNaN(mInput))
                .isEqualTo(mOutput);
    }    
}
