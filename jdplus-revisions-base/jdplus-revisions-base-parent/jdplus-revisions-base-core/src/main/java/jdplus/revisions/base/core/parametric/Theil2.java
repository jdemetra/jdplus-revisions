/*
 * Copyright 2020 National Bank of Belgium
 *
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved 
 * by the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 *
 * https://joinup.ec.europa.eu/software/page/eupl
 *
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and 
 * limitations under the Licence.
 */
package jdplus.revisions.base.core.parametric;

import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.data.DoubleSeqCursor;
import jdplus.toolkit.base.api.stats.StatException;
import lombok.NonNull;


/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
public class Theil2 {
    /**
     * !Temporary DUPLICATED from StatUtility in toolkit!
     * 
     * The second specification of Theil's inequality coefficient.  
     * This second specification avoid the problematic of near zero forecasts
     * of the first specification and has a clear interpretation.
     * 
     * <a href="https://www.economicsnetwork.ac.uk/showcase/cook_theil">Source</a>
     * 
     * @param a The first sequence
     * @param b The second sequence
     * @return The closer to 0, the more comparable the sequences.    
     */
    public double U2(@NonNull DoubleSeq a, @NonNull DoubleSeq b) {
        int n = a.length();
        
        if (b.length() != n) {
            throw new StatException("Non compatible data");
        }
        if(a.anyMatch(d -> d == 0)){
            return Double.NaN;
        }
        if(a.isEmpty() || b.isEmpty()){
            throw new StatException("a and b cannot be empty");
        }
        
        double nssq = 0, dssq = 0;
        DoubleSeqCursor atcur = a.cursor(), at1cur = a.cursor();
        DoubleSeqCursor bt1cur = b.cursor();
        at1cur.skip(1); bt1cur.skip(1);
       
        for (int i = 0; i < n-1; ++i) {
            double cat1 = at1cur.getAndNext(), cbt1 = bt1cur.getAndNext();
            double cat = atcur.getAndNext();
            double rn = cbt1 - cat1;
            rn /= cat;
            nssq += rn * rn;
            double rd = cat1 - cat;       
            rd /= cat;
            dssq += rd * rd;
        }
        if (dssq == 0) {
            return 0;
        }
        return Math.sqrt(nssq) / Math.sqrt(dssq);
    }
}


