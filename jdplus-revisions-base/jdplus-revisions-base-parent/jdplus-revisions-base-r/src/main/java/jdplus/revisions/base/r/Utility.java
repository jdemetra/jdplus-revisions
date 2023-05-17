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
package jdplus.revisions.base.r;

import java.time.LocalDate;
import jdplus.revisions.base.api.parametric.AutoCorrelationTests;
import jdplus.revisions.base.api.parametric.Bias;
import jdplus.revisions.base.api.parametric.Coefficient;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.data.DoubleSeqCursor;
import jdplus.revisions.base.api.parametric.OlsTests;
import jdplus.revisions.base.api.parametric.RegressionBasedAnalysis;
import jdplus.revisions.base.api.parametric.RevisionAnalysis;
import jdplus.revisions.base.api.parametric.SignalNoise;
import jdplus.revisions.base.api.parametric.UnitRoot;
import jdplus.revisions.base.core.parametric.AutoCorrelationTestsComputer;
import jdplus.revisions.base.core.parametric.BiasComputer;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.revisions.base.core.parametric.OlsTestsComputer;
import jdplus.revisions.base.core.parametric.SignalNoiseComputer;
import jdplus.toolkit.base.core.stats.StatUtility;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.revisions.base.core.treatment.PreTreatment;
import jdplus.revisions.base.core.parametric.Theil2;
import jdplus.revisions.base.core.parametric.UnitRootTestsComputer;
import jdplus.toolkit.base.api.data.DoublesMath;
import jdplus.toolkit.base.api.dstats.ContinuousDistribution;
import jdplus.toolkit.base.api.stats.StatisticalTest;
import jdplus.toolkit.base.api.stats.TestType;
import jdplus.toolkit.base.core.dstats.T;
import jdplus.toolkit.base.core.stats.tests.DickeyFuller;
import jdplus.toolkit.base.core.stats.tests.JohansenCointegration;
import jdplus.toolkit.base.core.stats.tests.TestsUtility;

/**
 *
 * @author PALATEJ
 */
@lombok.experimental.UtilityClass
public class Utility {

    /**
     * Theil coefficients computed on the columns of the vintages matrix
     *
     * @param vintages Vintages
     * @param gap Delay between the compared vintages (should be &ge 1)
     * @return
     */
    public double[] theil(Matrix vintages, int gap) {
        if (gap < 1) {
            throw new IllegalArgumentException("gap should be >= 1");
        }
        int n = vintages.getColumnsCount() - gap;
        if (n <= 0) {
            return null;
        }
        double[] u = new double[n];
        for (int i = 0; i < n; ++i) {
            DoubleSeq a = vintages.column(i + gap);
            DoubleSeq b = vintages.column(i);
            Matrix mc = PreTreatment.cleanNaN(a, b);
            u[i] = StatUtility.theilInequalityCoefficient(mc.column(0), mc.column(1));
        }
        return u;
    }

    /**
     * Theil2 coefficients computed on the columns of the vintages matrix
     *
     * @author clemasso
     *
     * @param vintages Vintages
     * @param gap Delay between the compared vintages (should be &ge 1)
     * @return
     */
    public double[] theil2(Matrix vintages, int gap) {
        if (gap < 1) {
            throw new IllegalArgumentException("gap should be >= 1");
        }
        int n = vintages.getColumnsCount() - gap;
        if (n <= 0) {
            return null;
        }
        double[] u = new double[n];
        for (int i = 0; i < n; ++i) {
            DoubleSeq a = vintages.column(i + gap);
            DoubleSeq b = vintages.column(i);
            Matrix mc = PreTreatment.cleanNaN(a, b);
            u[i] = Theil2.U2(mc.column(0), mc.column(1));
        }
        return u;
    }

    // apply it for other methods
    // use the modified code in R
    private final int OLS = 16, C = 3;

    /**
     * v(t)=a+b*v(t-gap)
     *
     * @param vintages Vintages
     * @param gap Delay between the compared vintages (should be &ge 1)
     * @return
     */
    public Matrix slopeAndDrift(Matrix vintages, int gap) {
        if (gap < 1) {
            throw new IllegalArgumentException("gap should be >= 1");
        }
        int n = vintages.getColumnsCount() - gap;
        if (n <= 0) {
            return null;
        }
        FastMatrix rslt = FastMatrix.make(n, OLS + 2 * C);

        for (int i = 0; i < n; ++i) {
            DoubleSeq y = vintages.column(i + gap);
            DoubleSeq x = vintages.column(i);
            Matrix yxCorr = PreTreatment.cleanNaN(y, x);
            DoubleSeqCursor.OnMutable cursor = rslt.row(i).cursor();
            OlsTests test = OlsTestsComputer.of(yxCorr.column(0), yxCorr.column(1));
            olsInformation(test, cursor);
            
            // Test beta1=1 instead of beta1=0
            double N = rslt.get(i, 0);
            int nx = 2;
            double slopeEst = rslt.get(i, 6);
            double slopeStdErr = rslt.get(i, 7);
            double t1 = (slopeEst-1)/slopeStdErr;
            T tdist = new T(N-nx);
            double pvalT1 = TestsUtility.pvalue(tdist, t1, TestType.TwoSided);
            rslt.set(i, 8, pvalT1);
        }
        
        return rslt;
    }

    private final int AC = 5;

    /**
     * v(t)=a+b*v(t-gap)
     *
     * @param vintages Vintages
     * @param nbg Number of lags in Breusch-Godfrey test
     * @param nlb Number of lag in Ljung-Box
     * @return
     */
    public Matrix autoCorrelation(Matrix vintages, int nbg, int nlb) {
        int n = vintages.getColumnsCount();
        FastMatrix rslt = FastMatrix.make(n * (n - 1) / 2, AC);

        for (int i = 0, k = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                try {
                    DoubleSeq y = vintages.column(i);
                    DoubleSeq x = vintages.column(j);
                    Matrix yxCorr = PreTreatment.cleanNaN(y, x);
                    DoubleSeqCursor.OnMutable cursor = rslt.row(k++).cursor();
                    AutoCorrelationTests test = AutoCorrelationTestsComputer.of(yxCorr.column(0), yxCorr.column(1), nbg, nlb);
                    acInformation(test, cursor);
                } catch (Exception err) {
                }
            }
        }
        return rslt;
    }

    private final int EG = 4;

    /**
     * v(t)=a+b*v(t-gap)
     *
     * @param vintages Vintages
     * @param adfk Number of lags in augmented dickey-fuller test
     * @return
     */
    public Matrix cointegration(Matrix vintages, int adfk) {
        int n = vintages.getColumnsCount();
        FastMatrix rslt = FastMatrix.make(n * (n - 1) / 2, EG);

        for (int i = 0, k = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                try {
                    DoubleSeq x = vintages.column(i);
                    DoubleSeq y = vintages.column(j);
                    Matrix xyCorr = PreTreatment.cleanNaN(x, y);

                    DoubleSeqCursor.OnMutable cursor = rslt.row(k++).cursor();
                    DickeyFuller df = DickeyFuller.engleGranger(xyCorr.column(1), xyCorr.column(0))
                            .numberOfLags(adfk).build();
                    if (df != null) {
                        cursor.setAndNext(df.getRho());
                        cursor.setAndNext(df.getSer());
                        cursor.setAndNext(df.getTest());
                        cursor.setAndNext(df.getPvalue());
                    }
                } catch (Exception err) {
                }
            }
        }
        return rslt;
    }

    private static final int JOHANSEN = 2;

    /**
     * v(t)=a+b*v(t-gap)
     *
     * @param vintages Vintages
     * @param lag Number of lags in augmented dickey-fuller test
     * @param model
     * @return
     */
    public Matrix vecm(Matrix vintages, int lag, String model) {
        int n = vintages.getColumnsCount();
        FastMatrix rslt = FastMatrix.make(n * (n - 1) / 2, JOHANSEN * lag);
        JohansenCointegration.ECDet ecdet = JohansenCointegration.ECDet.valueOf(model);
        JohansenCointegration computer = JohansenCointegration.builder()
                .errorCorrectionModel(ecdet)
                .lag(lag)
                .build();
        
        for (int i = 0, k = 0; i < n; ++i) { 
            for (int j = i + 1; j < n; ++j) {
                DoubleSeq vi = vintages.column(i);
                DoubleSeq vj = vintages.column(j);
                Matrix vijCorr = PreTreatment.cleanNaN(vi, vj);  
                FastMatrix M = FastMatrix.make(vijCorr.getRowsCount(), 2);
                M.column(0).copy(vijCorr.column(0));
                M.column(1).copy(vijCorr.column(1));
                try {
                    DoubleSeqCursor.OnMutable cursor = rslt.row(k++).cursor();
                    computer.process(M, null);
                    for (int l = lag - 1; l >= 0; --l) {
                        cursor.setAndNext(computer.traceTest(l));
                    }
                    for (int l = lag - 1; l >= 0; --l) {
                        cursor.setAndNext(computer.maxTest(l));
                    }
                } catch (Exception err) {
                }
            }
        }
        return rslt;
    }

    private final int UR = 4 * 4;

    /**
     * Computes unit roots tests. The tests are givenin the following order:
     * Dickey-Fuller Augmented Dickey-Fuller Dickey-Fuller with c and trend
     * Philips-Perron
     *
     * @param vintages
     * @param adfk Number of lags in augmented dickey-fuller test
     * @return
     */
    public Matrix unitroot(Matrix vintages, int adfk) {
        int n = vintages.getColumnsCount();
        FastMatrix rslt = FastMatrix.make(n, UR);

        for (int i = 0, k = 0; i < n; ++i) {
            try {
                DoubleSeqCursor.OnMutable cursor = rslt.row(k++).cursor();
                UnitRoot ur = UnitRootTestsComputer.of(PreTreatment.cleanNaN(vintages.column(i)), adfk);
                urInformation(ur, cursor);
            } catch (Exception err) {
            }
        }
        return rslt;
    }

    /**
     * rev(t)=a+b*v(t-gap)
     *
     * @param vintages Vintages
     * @param gap Delay between the compared vintages (should be &ge 1)
     * @return
     */
    public Matrix efficiencyModel1(Matrix vintages, int gap) {
        if (gap < 1) {
            throw new IllegalArgumentException("gap should be >= 1");
        }
        int n = vintages.getColumnsCount() - gap;
        if (n <= 0) {
            return null;
        }
        FastMatrix rslt = FastMatrix.make(n, OLS + 2 * C);

        for (int i = 0; i < n; ++i) {
            
            DoubleSeq x = vintages.column(i);
            DoubleSeq y = DoublesMath.subtract(vintages.column(i + gap), x);
            Matrix yxCorr = PreTreatment.cleanNaN(y, x);
            
            DoubleSeqCursor.OnMutable cursor = rslt.row(i).cursor();
            OlsTests test = OlsTestsComputer.of(yxCorr.column(0), yxCorr.column(1));
            olsInformation(test, cursor);
        }
        return rslt;
    }

    /**
     * rev(t)=a+b*rev(t-1)
     *
     * @param vintages Vintages
     * @param gap Delay between the vintages used to compute the revisions
     * (should be &ge 1)
     * @return
     */
    public Matrix efficiencyModel2(Matrix vintages, int gap) {
        int n = vintages.getColumnsCount() - gap - 1;
        FastMatrix rslt = FastMatrix.make(n, OLS + 2 * C);
        for (int i = 0; i < n; ++i) {
            DoubleSeqCursor.OnMutable cursor = rslt.row(i).cursor();
            try {
                DoubleSeq y = DoublesMath.subtract(vintages.column(i + gap + 1), vintages.column(i + 1));
                DoubleSeq x = DoublesMath.subtract(vintages.column(i + gap), vintages.column(i));
                Matrix yxCorr = PreTreatment.cleanNaN(y, x);
                
                OlsTests test = OlsTestsComputer.of(yxCorr.column(0), yxCorr.column(1));
                olsInformation(test, cursor);
            } catch (Exception err) {
            }
        }
        return rslt;
    }

    /**
     * rev(t)=a+b(1)*rev(t-1)+b(2)*rev(t-2)+...+b(nrevs)*rev(t-nrevs)
     *
     * @param revs
     * @param nrevs
     * @return
     */
    public Matrix orthogonallyModel1(Matrix revs, int nrevs) {
        int nr = revs.getRowsCount();
        int nc = revs.getColumnsCount();
        if (nrevs >= nc) {
            return null;
        }
        FastMatrix rslt = FastMatrix.make(nc - nrevs, OLS + C * (1 + nrevs));
        DoubleSeq[] x = new DoubleSeq[nrevs];
        for (int i = nrevs; i < nc; ++i) {
            double[] yx = new double[nr * (nrevs + 1)];
            double[] y = revs.column(i).toArray();
            System.arraycopy(y, 0, yx, 0, y.length);

            for (int j = 0; j < nrevs; ++j) {
                x[j] = revs.column(i - j - 1);
                System.arraycopy(x[j].toArray(), 0, yx, (j + 1) * nr, x[j].toArray().length);
            }

            Matrix yxCorr = PreTreatment.cleanNaN(Matrix.of(yx, nr, nrevs + 1));

            DoubleSeq yc = yxCorr.column(0);
            DoubleSeq[] xc = new DoubleSeq[nrevs];
            for (int k = 0; k < nrevs; ++k) {
                xc[k] = yxCorr.column(k + 1);
            }

            DoubleSeqCursor.OnMutable cursor = rslt.row(i - nrevs).cursor();
            try {
                OlsTests test = OlsTestsComputer.of(yc, xc);
                olsInformation(test, cursor);
            } catch (Exception err) {
            }
        }
        return rslt;
    }
    
    public Matrix orthogonallyModel2(Matrix revs, int k) {
        int n = revs.getColumnsCount();
        if (k >= n || k < 1) {
            return null;
        }
        FastMatrix rslt = FastMatrix.make(n - k, OLS + C * 2);

        for (int i = 0, j = 0; i < n; ++i) {
            if (i > k - 1) {

                DoubleSeq y = revs.column(i);
                DoubleSeq x = revs.column(i - k);
                Matrix yxCorr = PreTreatment.cleanNaN(y, x);

                DoubleSeqCursor.OnMutable cursor = rslt.row(j++).cursor();
                try {
                    OlsTests test = OlsTestsComputer.of(yxCorr.column(0), yxCorr.column(1));
                    olsInformation(test, cursor);
                } catch (Exception err) {
                }
            }
        }
        return rslt;
    }

    public double theil(RegressionBasedAnalysis<LocalDate> analysis, int k) {
        if (k > analysis.getRevisions().size()) {
            return Double.NaN;
        }
        return analysis.getRevisions().get(k - 1).getTheilCoefficient();
    }

    private static final int BIAS = 9;

    /**
     * Bias computed on a matrix of revisions (each column corresponds to a
     * revision)
     *
     * @param revs The revisions
     * @return
     */
    public Matrix bias(Matrix revs) {
        int n = revs.getColumnsCount();
        FastMatrix rslt = FastMatrix.make(n, BIAS);

        for (int i = 0; i < n; ++i) {
            DoubleSeq cur = revs.column(i);
            DoubleSeqCursor.OnMutable cursor = rslt.row(i).cursor();
            Bias bias = BiasComputer.of(cur);
            biasInformation(bias, cursor);
        }
        return rslt;
    }

    private final int SN = 6;

    public Matrix signalNoise(Matrix vintages, int gap) {
        if (gap < 1) {
            throw new IllegalArgumentException("gap should be >= 1");
        }
        int n = vintages.getColumnsCount() - gap;
        if (n <= 0) {
            return null;
        }
        FastMatrix rslt = FastMatrix.make(n, SN);

        for (int i = 0; i < n; ++i) {
            DoubleSeq L = vintages.column(i + gap);
            DoubleSeq P = vintages.column(i);
            Matrix LPCorr = PreTreatment.cleanNaN(L, P);

            DoubleSeqCursor.OnMutable cursor = rslt.row(i).cursor();
            SignalNoise test = SignalNoiseComputer.of(LPCorr.column(1), LPCorr.column(0));
            signalNoiseInformation(test, cursor);
        }
        return rslt;
    }

    public void olsInformation(OlsTests reg, DoubleSeqCursor.OnMutable cursor) {
        if (reg == null) {
            return;
        }
        Coefficient[] c = reg.getCoefficients();
        StatisticalTest jb = reg.getDiagnostics().getJarqueBera();
        StatisticalTest bp = reg.getDiagnostics().getBreuschPagan();
        StatisticalTest w = reg.getDiagnostics().getWhite();
        StatisticalTest arch = reg.getDiagnostics().getArch();
        cursor.setAndNext(reg.getN());
        cursor.setAndNext(reg.getR2());
        cursor.setAndNext(reg.getF());
        for (int i = 0; i < c.length; ++i) {
            cursor.setAndNext(c[i].getEstimate());
            cursor.setAndNext(c[i].getStdev());
            cursor.setAndNext(c[i].getPvalue());
        }
        cursor.setAndNext(reg.getDiagnostics().getSkewness());
        cursor.setAndNext(reg.getDiagnostics().getKurtosis());
        cursor.setAndNext(jb.getValue());
        cursor.setAndNext(jb.getPvalue());
        cursor.setAndNext(reg.getDiagnostics().getBpr2());
        cursor.setAndNext(bp.getValue());
        cursor.setAndNext(bp.getPvalue());
        cursor.setAndNext(reg.getDiagnostics().getWr2());
        cursor.setAndNext(w.getValue());
        cursor.setAndNext(w.getPvalue());
        cursor.setAndNext(reg.getDiagnostics().getArchr2());
        cursor.setAndNext(arch.getValue());
        cursor.setAndNext(arch.getPvalue());
    }

    public void acInformation(AutoCorrelationTests ac, DoubleSeqCursor.OnMutable cursor) {
        if (ac == null) {
            return;
        }
        StatisticalTest bg = ac.getBreuschGodfrey();
        StatisticalTest lb = ac.getLjungBox();
        cursor.setAndNext(ac.getBgr2());
        cursor.setAndNext(bg.getValue());
        cursor.setAndNext(bg.getPvalue());
        cursor.setAndNext(lb.getValue());
        cursor.setAndNext(lb.getPvalue());
    }

    public double[] biasInformation(RegressionBasedAnalysis<LocalDate> analysis, int k) {
        if (k > analysis.getRevisions().size()) {
            return null;
        }
        RevisionAnalysis<LocalDate> cur = analysis.getRevisions().get(k - 1);
        if (cur == null) {
            return null;
        }
        Bias bias = cur.getBias();
        if (bias == null) {
            return null;
        }
        return new double[]{
            bias.getN(),
            bias.getMu(),
            bias.getSigma(),
            bias.getT(),
            bias.getTPvalue(),
            bias.getAr(),
            bias.getAdjustedSigma(),
            bias.getAdjustedT(),
            bias.getAdjustedTPvalue()};
    }

    public void biasInformation(Bias bias, DoubleSeqCursor.OnMutable cursor) {
        if (bias == null) {
            return;
        }
        cursor.setAndNext(bias.getN());
        cursor.setAndNext(bias.getMu());
        cursor.setAndNext(bias.getSigma());
        cursor.setAndNext(bias.getT());
        cursor.setAndNext(bias.getTPvalue());
        cursor.setAndNext(bias.getAr());
        cursor.setAndNext(bias.getAdjustedSigma());
        cursor.setAndNext(bias.getAdjustedT());
        cursor.setAndNext(bias.getAdjustedTPvalue());
    }

    private static void urInformation(UnitRoot ur, DoubleSeqCursor.OnMutable cursor) {
        if (ur == null) {
            return;
        }
        UnitRoot.Test t = ur.getDickeyFuller();
        cursor.setAndNext(t.getValue());
        cursor.setAndNext(t.getStdev());
        cursor.setAndNext(t.getStatistic());
        cursor.setAndNext(t.getPvalue());
        t = ur.getAugmentedDickeyFuller();
        cursor.setAndNext(t.getValue());
        cursor.setAndNext(t.getStdev());
        cursor.setAndNext(t.getStatistic());
        cursor.setAndNext(t.getPvalue());
        t = ur.getDickeyFullerWithTrendAndIntercept();
        cursor.setAndNext(t.getValue());
        cursor.setAndNext(t.getStdev());
        cursor.setAndNext(t.getStatistic());
        cursor.setAndNext(t.getPvalue());
        t = ur.getPhilipsPerron();
        cursor.setAndNext(t.getValue());
        cursor.setAndNext(t.getStdev());
        cursor.setAndNext(t.getStatistic());
        cursor.setAndNext(t.getPvalue());
    }

    private static void signalNoiseInformation(SignalNoise sn, DoubleSeqCursor.OnMutable cursor) {
        if (sn == null) {
            return;
        }
        cursor.setAndNext(sn.getNewsR2());
        cursor.setAndNext(sn.getNewsF());
        cursor.setAndNext(sn.getNewsPvalue());
        cursor.setAndNext(sn.getNoiseR2());
        cursor.setAndNext(sn.getNoiseF());
        cursor.setAndNext(sn.getNoisePvalue());
    }
    
    
    
}
