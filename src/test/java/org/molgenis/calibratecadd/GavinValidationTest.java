package org.molgenis.calibratecadd;

import org.testng.annotations.Test;

/**
 * Created by joeri on 9/1/16.
 */
public class GavinValidationTest {

    @Test
    public void testPredictionTool() throws Exception {
        new AbstractValidationTest("GAVIN").testPredictionTool();
    }

}
