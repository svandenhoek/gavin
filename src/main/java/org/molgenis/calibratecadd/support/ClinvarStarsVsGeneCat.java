package org.molgenis.calibratecadd.support;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.molgenis.data.annotation.entity.impl.gavin.GavinEntry;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;

/**
 * meanreviewstatus per gene calculator
 * Created by joeri on 1/17/17.
 */
public class ClinvarStarsVsGeneCat {

    public ClinvarStarsVsGeneCat(String[] args) throws Exception {
        HashMap<String, ArrayList<String>> genes = new HashMap<String, ArrayList<String>>();
        // GAVIN_calibrations_r0.3.tsv
        GavinUtils gu = new GavinUtils(new File(args[0]));

        for(String gene : gu.geneToEntry.keySet())
        {
            GavinEntry ge = gu.geneToEntry.get(gene);

         //   if(ge.category == GavinEntry.Category.C1 || ge.category == GavinEntry.Category.C2 || ge.category == GavinEntry.Category.C4)
         //   {
                genes.put(gene, new ArrayList<>());
        //    }
        }


        // variant_summary.txt
        Scanner vs = new Scanner(new File(args[1]));

        //skip header
        String line = vs.nextLine();
        System.out.println("header: " + line);

        while (vs.hasNextLine()) {
            line = vs.nextLine();
            String[] lineSplit = line.split("\t", -1);
            String genomeBuild = lineSplit[16];

            if (!genomeBuild.equals("GRCh37") && !genomeBuild.equals("GRCh38") && !genomeBuild.equals("NCBI36") & !genomeBuild.equals("")) {
                throw new Exception("bad genome build: " + genomeBuild);
            }
            // needs to be GRCh37
            if (!genomeBuild.equals("GRCh37")) {
                continue;
            }

            // see: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/README
            String GeneSymbol = lineSplit[4];
            String ReviewStatus = lineSplit[24];
            String NumberSubmitters = lineSplit[25];
            String Guidelines = lineSplit[26];
            String SubmitterCategories = lineSplit[29];


            if(genes.keySet().contains(GeneSymbol))
            {
                genes.get(GeneSymbol).add(ReviewStatusToNr(ReviewStatus) + "_" + Short.parseShort(NumberSubmitters.equals("")?"1":NumberSubmitters) + "_" + Short.parseShort(SubmitterCategories));
            }
        }





        System.out.println("Gene" + "\t" + "Category" + "\t" + "Pvalue" + "\t" +"NrOfVariants" + "\t" + "ReviewStatusMean" + "\t" + "NumberSubmittersMean" + "\t" + "SubmitterCategoriesMean");

        for(String gene : genes.keySet())
        {
            if(genes.get(gene).size() == 0)
            {
                continue;
            }

            double[] ReviewStatusArr = new double[genes.get(gene).size()];
            double[] NumberSubmittersArr = new double[genes.get(gene).size()];
            double[] SubmitterCategoriesArr = new double[genes.get(gene).size()];


            for(int i = 0; i < genes.get(gene).size(); i++)
            {
                String[] split = genes.get(gene).get(i).split("_");
                ReviewStatusArr[i] = Double.parseDouble(split[0]);
                NumberSubmittersArr[i] = Double.parseDouble(split[1]);
                SubmitterCategoriesArr[i] = Double.parseDouble(split[2]);

            }

            Mean mean = new Mean();
            double ReviewStatusMean = mean.evaluate(ReviewStatusArr);
            double NumberSubmittersMean = mean.evaluate(NumberSubmittersArr);
            double SubmitterCategoriesMean = mean.evaluate(SubmitterCategoriesArr);

            System.out.println(gene + "\t" + gu.geneToEntry.get(gene).category + "\t" + gu.geneToEntry.get(gene).UTestPvalue + "\t" + genes.get(gene).size() + "\t" + ReviewStatusMean + "\t" + NumberSubmittersMean + "\t" + SubmitterCategoriesMean);
        }





    }

    /**
     * see: https://www.ncbi.nlm.nih.gov/clinvar/docs/variation_report/#interpretation
     * @param ReviewStatus
     * @return
     * @throws Exception
     */
    public Short ReviewStatusToNr(String ReviewStatus) throws Exception {
        if(ReviewStatus.equals("no assertion criteria provided"))
        {
            return 0;
        }
        if(ReviewStatus.equals("no assertion for the individual variant"))
        {
            return 0;
        }
        if(ReviewStatus.equals("no assertion provided"))
        {
            return 0;
        }
        if(ReviewStatus.equals("criteria provided, single submitter"))
        {
            return 1;
        }
        if(ReviewStatus.equals("criteria provided, conflicting interpretations"))
        {
            return 1;
        }
        if(ReviewStatus.equals("criteria provided, multiple submitters, no conflicts"))
        {
            return 2;
        }
        if(ReviewStatus.equals("reviewed by expert panel"))
        {
            return 3;
        }
        if(ReviewStatus.equals("practice guideline"))
        {
            return 4;
        }
        throw new Exception("unknown ReviewStatus: " + ReviewStatus);
    }

    public static void main (String[] args) throws Exception {

        new ClinvarStarsVsGeneCat(args);

    }
}
