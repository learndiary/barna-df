package barna.geneid;

import barna.commons.log.Log;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 10:41 AM
 * To change this template use File | Settings | File Templates.
 */
public class GeneID {

    /**
     * Matrix to translate characters to numbers. borrowed from jwf
     */
    public static final int TRANS[] = {
            /* Control characters */
            4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,
            /* Punctuation and digits */
            4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,
            /* Capitals */     /*  A=0; C=1; G=2; T=3; other = 4  */
            4,0,4,1,4,4,4,2,   /* @,A-G: A,C and G found */
            4,4,4,4,4,4,4,4,   /* H-O   */
            4,4,4,4,3,3,4,4,   /* P-W: T and U found */
            4,4,4,4,4,4,4,4,   /* X-Z,etc */
            /* Lower case */
            4,0,4,1,4,4,4,2,   /*  @,A-G  */
            4,4,4,4,4,4,4,4,   /*   H-O   */
            4,4,4,4,3,3,4,4,   /*   P-W   */
            4,4,4,4,4,4,4,4    /* X-Z,etc */
    };

    /**
     * Length of every processed fragment
     */
    public static final int LENGTHSi= 220000;

    /**
     * One signal per L / RSITES bp
     */
    public static final int RSITES= 1;

    /**
     * One exon per L / REXONS bp
     */
    public static final int REXONS= 2;    // was 3

    /**
     * Estimated amount of backup signals
     */
    public static final int RBSITES= 75;

    /**
     * Estimated amount of backup exons
     */
    public static final int RBEXONS= 125; // was 125


    /* Basic values (in addition to ratios)     */

    public static final int BASEVALUESITES_SHORT= 100000;
    public static final int BASEVALUEEXONS_SHORT= 6000;
    public static final int BASEVALUESITES_LARGE= 600000;
    public static final int BASEVALUEEXONS_LARGE= 600000;


    /* Generic maximum values: sites, exons and backup elements */

    /**
     * From the defined values RSITES (geneid.h) and length of sequence,
     * an estimation for the amount of predicted signals and exons (any type) is
     * computed in order to ask for enough memory to allocate them
     */
    long NUMSITES;

    /**
    * From the defined values REXONS (geneid.h) and length of sequence,
    * an estimation for the amount of predicted signals and exons (any type) is
    * computed in order to ask for enough memory to allocate them
    */
    long NUMEXONS;

    /**
     * From the defined values RBSITES (geneid.h) and length of sequence,
     * an estimation for the amount of signals and exons (any type), necessary
     * to restore the prediction between 2 splits, is computed in order to ask
     * for enough memory to allocate them
     */
    long MAXBACKUPSITES;

    /**
     * From the defined values RBEXONS (geneid.h) and length of sequence,
     * an estimation for the amount of signals and exons (any type), necessary
     * to restore the prediction between 2 splits, is computed in order to ask
     * for enough memory to allocate them
     */
    long MAXBACKUPEXONS;

    long NUMU12SITES,NUMU12EXONS,NUMU12U12EXONS;



    /**
     * Sets the generic maximum values (sites, exons and backup elements)
     * based on the observed sequence (chromosome) length
     * @param L
     */
    void setRatios(long L) {

        // TODO move to constructor

        /* L is the estimated length of input DNA sequence */
        /* LENGTHSi is the length splitting */
        if (L < LENGTHSi)
        {
            /* Short sequences processed as a whole: only one split */
            NUMSITES = L / RSITES + BASEVALUESITES_SHORT;
            NUMEXONS = L / REXONS + BASEVALUEEXONS_SHORT;

            /* There is no need to divide the sequence */
            MAXBACKUPSITES = 0;
            MAXBACKUPEXONS = 0;
        }
        else
        {
            /* Long sequences must be divided into several fragments */
            NUMSITES = LENGTHSi / RSITES;
            NUMEXONS = LENGTHSi / REXONS;

            /* Information inter-split predictions must be saved */
            MAXBACKUPSITES = (L / RBSITES) + BASEVALUESITES_LARGE;
            MAXBACKUPEXONS = (L / RBEXONS) + BASEVALUEEXONS_LARGE;
        }

    }


    /**
     * Translation from string into integer: for oligonucleotides with this length
     *
     * @param s sequence
     * @param ls length of sequence
     * @param cardinal length of the alphabet
     * @return
     */
    static int OligoToInt(String s, int ls, int cardinal) {

        int index = 0;
        int weight = 1;

        // index = 5*TRANS[(int)(*(s + i -1))] + TRANS[(int)(*(s + i))];
        for (int i= ls- 1; i>= 0; --i) {
            index += weight* TRANS[(int) s.charAt(i)];
            weight *= cardinal;
        }

        return index;
    }


    /* Search for acceptor splice sites, using additional profiles */
    long  BuildAcceptors(String s,
                         short klass,
                         String type,
                         String subtype,
                         Profile p,
                         Profile ppt,
                         Profile bp,
                         Site[] st,
                         int l1,
                         int l2,
                         int ns,
                         int nsites,
                         int Strand,
                         PackExternalInformation external) {
        int i,j;
        String sOriginal;
        double score;
        double scoreBP;
        double scorePPT;
        double scoreAcc;
        /*   long ns,is; */
        int is;
        int left,right;
        int index;
        float cutoff;
        /* Final number of predicted signals (that type) */
        /*   ns = 0; */

        /* Back-up the origin of the sequence */
        sOriginal = s;

        cutoff = p.cutoff; /* For U2 acceptors, we currently use cutoff given in param file*/


        /* 1. Searching sites between beginning of the sequence and p.offset */
        if (l1== 0)
        {
            for (is = 0; is < (p.offset)  && (ns< NUMSITES); is++)
            {
                score= 0f;
                /* Applying part of the profile */
                for (i=p.offset- is, j=0; i < p.dimension; i++,j++)
                {
                    /* i is the position inside the region */
                    index = OligoToInt(s.substring(9), p.order+ 1,5);   // s+ j

                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }

                scorePPT = 0f;
                scoreBP = 0f;
                scoreAcc = 0f;
                if (score >= cutoff){
                    /* Using additional profiles */
                    if (GParam.PPT!= 0)
                        scorePPT = computeU2PPTProfile(sOriginal,p.offset-is,l2,ppt, st[ns]);

                    if (GParam.BP!= 0)
                        scoreBP = computeU2BranchProfile(sOriginal,p.offset-is,l2,bp, st[ns]);

                    /* For the time being, we will not use the BP or PPT scores */
                    /* if (scoreBP > 0){score = score + scoreBP;} */ /* + scorePPT */
                    scoreAcc = score;
                    score = score + scoreBP;
                    score = p.afactor + (p.bfactor * score);
                    if(GeneIDsettings.UTR!= 0){
                        score = score + peakEdgeScore(is + p.order, Strand, external, l1, l2, 6);
                    }
                    /* Acceptor core is used as a global cutoff */

                    if (score >= p.cutoff)
                    {
                        st[ns].Position = is + p.order;
                        st[ns].ScoreBP = scoreBP;
                        st[ns].ScorePPT = scorePPT;
                        st[ns].ScoreAccProfile = scoreAcc;
                        st[ns].Score = score;
                        st[ns].sclass= klass;
                        st[ns].type= type;
                        st[ns].subtype= subtype;
                        ns++;
                    }
                }
            }
        }

        /* 2. Normal processing: predicting using the whole profile */
        /* left and right are the true boundaries of prediction */
        left  = Math.max(0+ p.order, l1 - p.offset); // MAX()
        right = l2 - p.offset;
        s += left;
        is = 0;
        /* Case A: Using Markov chain with order 0: PWM */
        if (p.order == 0)
        {
            /* discovering splice sites with current profile */
            while (s.length()< (p.dimension- 1) && (is < right- left + 1) && (ns< NUMSITES)) // *(s+ p.dimension- 1) &&..
            {
                /* is = 0..right */
                score=0f;
                for (i= 0; i< p.dimension; i++)
                {
                    /* i is the position inside the region */
                    index = TRANS[(int) s.charAt(i)];    // (*(s + i))
                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }

                scorePPT = 0f;
                scoreBP = 0f;
                scoreAcc = 0f;
                if (score >= cutoff){
                    /* Using additional profiles */
                    if (GParam.PPT!= 0)
                        scorePPT = computeU2PPTProfile(sOriginal,left + is + p.offset,l2,ppt, st[ns]);
                    if (GParam.BP!= 0)
                        scoreBP = computeU2BranchProfile(sOriginal,left + is + p.offset,l2,bp, st[ns]);

                    /* if (scoreBP > 0){score = score + scoreBP;} */ /* + scorePPT */
                    scoreAcc = score;
                    score = score + scoreBP;
                    score = p.afactor + (p.bfactor * score);
                    if(GeneIDsettings.UTR!= 0){
                        score = score + peakEdgeScore(left + is + p.offset, Strand, external, l1, l2, 6);
                    }
                    if (score >= p.cutoff)
                    {
                        st[ns].Position = left + is + p.offset;
                        st[ns].ScoreBP = scoreBP;
                        st[ns].ScorePPT = scorePPT;
                        st[ns].sclass= klass;
                        st[ns].ScoreAccProfile = scoreAcc;
                        st[ns].Score = score;
                        st[ns].type= type;
                        st[ns].subtype= subtype;
                        ns++;
                    }
                }
                is++;
                s= s.substring(1);  // s++;
            }
        }
        /* case B: Using Markov chain with order 1: dinucleotides */
        else if (p.order == 1)
        {
            /* discovering splice sites with current profile */
            while (s.length()< (p.dimension- 1) && (is < right- left + 1) && (ns<NUMSITES))  // *(s+p.dimension -1) && ..
            {
                /* is = 0..right */
                score=0f;
                for (i=0;i<p.dimension;i++)
                {
                    /* i is the position inside the region */
                    index = 5* TRANS[(int) s.charAt(i- 1)] + TRANS[(int) s.charAt(i)];
                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }

                scorePPT = 0f;
                scoreBP = 0f;
                scoreAcc = 0f;
                if (score >= cutoff){
                    /* Using additional profiles */
                    if (GParam.PPT!= 0)
                        scorePPT = computeU2PPTProfile(sOriginal,left + is + p.offset,l2,ppt, st[ns]);

                    if (GParam.BP!= 0)
                        scoreBP = computeU2BranchProfile(sOriginal,left + is + p.offset,l2,bp, st[ns]);

                    /* if (scoreBP > 0){score = score + scoreBP;} */ /* + scorePPT */
                    scoreAcc = score;
                    score = score + scoreBP;
                    score = p.afactor + (p.bfactor * score);
                    if(GeneIDsettings.UTR!= 0){
                        score = score + peakEdgeScore(left + is + p.offset, Strand, external, l1, l2, 6);
                    }

                    if (score >= p.cutoff)
                    {
                        st[ns].Position = left + is + p.offset;
                        st[ns].ScoreBP = scoreBP;
                        st[ns].ScorePPT = scorePPT;
                        st[ns].ScoreAccProfile = scoreAcc;
                        st[ns].Score = score;
                        st[ns].sclass= klass;
                        st[ns].type= type;
                        st[ns].subtype= subtype;
                        ns++;
                    }
                }
                is++;
                s= s.substring(1);  // s++
            }
        }
        /* case C: Using Markov chain with order > 1 */
        else
        {
            /* discovering splice sites with current profile */
            while (s.length()< (p.dimension- 1) && (is < right- left + 1) && (ns<NUMSITES)) // *(s+p.dimension -1) && ..
            {
                /* is = 0..right */
                score= 0f;
                for (i=0;i<p.dimension;i++)
                {
                    /* i is the position inside the region */
                    /* 5 is used because there are A,C,G,T and N */
                    index = OligoToInt(s.substring(i - p.order), p.order+1,5);

                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }

                scorePPT = 0f;
                scoreBP = 0f;
                scoreAcc = 0f;
                if (score >= cutoff){
                    /* Using additional profiles */
                    if (GParam.PPT!= 0)
                        scorePPT = computeU2PPTProfile(sOriginal,left + is + p.offset,l2,ppt, st[ns]);

                    if (GParam.BP!= 0)
                        scoreBP = computeU2BranchProfile(sOriginal,left + is + p.offset,l2,bp, st[ns]);

                    /* if (scoreBP > 0){score = score + scoreBP;} */ /* + scorePPT */
                    scoreAcc = score;
                    score = score + scoreBP;
                    score = p.afactor + (p.bfactor * score);
                    if(GeneIDsettings.UTR!= 0){
                        score = score + peakEdgeScore(left + is + p.offset, Strand, external, l1, l2, 6);
                    }

                    if (score >= p.cutoff)
                    {
                        st[ns].Position = left + is + p.offset;
                        st[ns].ScoreBP = scoreBP;
                        st[ns].ScorePPT = scorePPT;
                        st[ns].ScoreAccProfile = scoreAcc;
                        st[ns].Score = score;
                        st[ns].sclass= klass;
                        st[ns].type= type;
                        st[ns].subtype= subtype;
                        ns++;
                    }
                }
                is++;
                s= s.substring(1);  // s++
            }
        }
        if (ns >= nsites)
            Log.error("Too many predicted sites: decrease RSITES parameter");

        return(ns);
    }

    public static float scoreAcceptor(String s, Profile profileAcc, Profile profilePPT, Profile profileBP) {

        float score= scoreSite(s, profileAcc);

//        float scorePPT = 0f;
//        float scoreBP = 0f;
//        float scoreAcc = 0f;
/*        if (score >= profileAcc.cutoff){
            // Using additional profiles
            if (GParam.PPT!= 0)
                scorePPT = computeU2PPTProfile(sOriginal,profileAcc.offset-is,l2,ppt, st[ns]);

            if (GParam.BP!= 0)
                scoreBP = computeU2BranchProfile(sOriginal,profileAcc.offset-is,l2,bp, st[ns]);

            // For the time being, we will not use the BP or PPT scores
            // if (scoreBP > 0){score = score + scoreBP;} */ /* + scorePPT
            scoreAcc = score;
            score = score + scoreBP;
            score = profileAcc.afactor + (profileAcc.bfactor * score);
            if(GeneIDsettings.UTR!= 0){
                score = score + peakEdgeScore(is + profileAcc.order, Strand, external, l1, l2, 6);
            }
            // Acceptor core is used as a global cutoff

            if (score >= profileAcc.cutoff)
            {
                st[ns].Position = is + profileAcc.order;
                st[ns].ScoreBP = scoreBP;
                st[ns].ScorePPT = scorePPT;
                st[ns].ScoreAccProfile = scoreAcc;
                st[ns].Score = score;
                st[ns].sclass= klass;
                st[ns].type= type;
                st[ns].subtype= subtype;
                ns++;
            }
        }
*/
        return score;
    }

    /**
     * Returns the sum of log-likelihood scores for a certain sequence
     * matching the given site profile.
     * @param s sequence that is evaluated
     * @param p site profile for evaluation
     * @return the log-likelihood score
     */
    public static float scoreSite(String s, Profile p) {
        float score= 0f;
        if (s.length()!= p.dimension ) //+ p.order)
            throw new RuntimeException("Site sequence "+ s.length()+ "nt, expected "+ p.dimension+ "nt!");

        /* Applying part of the profile */
        //System.err.println(p.offset+",dim="+p.dimension+",slen="+s.length());
        //for (int i= 0; i < p.dimension; i++) {
        for (int i= 0; i < (p.dimension- p.order- 1); i++) {

            /* i is the position inside the region */
            int index = OligoToInt(s.substring(i), p.order+ 1, 5);

            float inc= 0;
            if (index >= p.dimensionTrans)
                inc= -GeneIDconstants.INFI;
            else
                //inc= p.transitionValues[i][index]; // reverted in lncc_fix, this line was probably wrong
                inc= p.transitionValues[i+ 1][index];
            score+= inc;
        }
        score = p.afactor + (p.bfactor * score);

        return score;
    }


    public static float scoreDonor(String s, Profile p) {

        // TODO this is a simplified version of the BuildDonors.c

        float score= 0f;
        // prefix_donor = (offset+ order)
        // prefix_acceptor = (dimension- offset- 2)+ order
        // suffix_donor = ..
        // suffix_acceptor = ..


        /* Applying part of the profile */
        for (int i=p.offset, j=0; i < p.dimension; i++, j++) {

            /* i is the position inside the region */
            String tuple= s.substring(j);   // works for GCAG
            int index = OligoToInt(tuple, p.order+ 1, 5);

            float incr= 0;
            if (index >= p.dimensionTrans)
                incr= -GeneIDconstants.INFI;
            else
                // incr= p.transitionValues[j][index]; // reverted in lncc_fix, this line was probably wrong
                incr= p.transitionValues[i][index];
            score+= incr;
        }
        score = p.afactor + (p.bfactor * score);

        return score;
    }



    /**
     *
     * @param s sequence to be scored
     * @param klass copied to each site found, TODO REMOVE
     * @param type copied to each site found, TODO REMOVE
     * @param subtype copied to each site found, TODO REMOVE
     * @param p donor site profile
     * @param st array to store objects that describe the site predictions
     * @param l1 flag, switches with <code>0</code> to "normal" prediction mode (using the whole profile),
     *           otherwise sites are searched between the beginning of the sequence and <code>p.offset</code>
     * @param l2 sequence clipping, marks the end of prediction in the string provided as <code>s</code>
     * @param ns number of sites that have already been predicted
     * @param nsites maximum number of sites to predict
     * @param strand directionality of prediction
     * @param external
     * @return
     */
    long buildDonors(String s,
                      short klass,
                      String type,
                      String subtype,
                      Profile p,
                      Site[] st,
                      int l1,
                      int l2,
                      int ns,
                      long nsites,
                      int strand,
                      PackExternalInformation external) {

        float score;
        int left,right;
        int index;

        /* 1. Searching sites between beginning of the sequence and p.offset */
        if (l1!= 0) {

            for (int is = 0; is < p.offset && (ns<nsites); is++)
            {
                score=0f;
                /* Applying part of the profile */
                for (int i= p.offset- is, j=0; i < p.dimension; i++,j++)
                {
                    /* i is the position inside the region */
                    index = OligoToInt(s.substring(j), p.order+ 1,5);

                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }
                score = p.afactor + (p.bfactor * score);
                if(GeneIDsettings.UTR!= 0){
                    score = score - peakEdgeScore(is + p.order, strand, external, l1, l2, 6);
                }
                if (score >= p.cutoff)
                {
                    st[ns].Position = is + p.order;
                    st[ns].Score= score;
                    st[ns].sclass= klass;
                    st[ns].type= type;
                    st[ns].subtype= subtype;
                    ns++;
                }
            }
        }

        /* 2. Normal processing: predicting using the whole profile */
        /* left and right are the true boundaries of prediction */
        left  = Math.max(0+ p.order, l1- p.offset); // MAX()
        right = l2 - p.offset;
        s += left;
        int is = 0;
        /* Case A: Using Markov chain with order 0: PWM */
        if (p.order == 0)
        {
            /* discovering splice sites with current profile */
            while (s.length()< (p.dimension- 1) && (is < right- left + 1) && (ns<nsites))   // *(s+p.dimension-1) &&...
            {
                /* is = 0..right */
                score=0f;
                for (int i=0; i< p.dimension; i++)
                {
                    /* i is the position inside the region */
                    index = TRANS[(int) s.charAt(i)];   // *(s + i))
                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }
                score = p.afactor + (p.bfactor * score);
                if(GeneIDsettings.UTR!= 0){
                    score = score - peakEdgeScore(left + is + p.offset, strand, external, l1, l2, 6);
                }
                if (score >= p.cutoff)
                {
                    st[ns].Position = left + is + p.offset;
                    st[ns].Score=score;
                    st[ns].sclass= klass;
                    st[ns].type= type;
                    st[ns].subtype= subtype;
                    ns++;
                }
                is++;
                s= s.substring(1);  // s++
            }
        }
        /* case B: Using Markov chain with order 1: dinucleotides */
        else if (p.order == 1)
        {
            /* discovering splice sites with current profile */
            while (s.length()< (p.dimension- 1) && (is < right- left + 1) && (ns<nsites))   // *(s+p.dimension-1) &&
            {
                /* is = 0..right */
                score=0f;
                for (int i=0;i<p.dimension;i++)
                {
                    /* i is the position inside the region */
                    index = 5*TRANS[(int) s.charAt(i- 1)] + TRANS[(int) s.charAt(i)];  // (*(s + i -1))  .. (*(s + i))
                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }
                score = p.afactor + (p.bfactor * score);
                if(GeneIDsettings.UTR!= 0){
                    score = score - peakEdgeScore(left + is + p.offset, strand, external, l1, l2, 6);
                }
                if (score >= p.cutoff)
                {
                    st[ns].Position = left + is + p.offset;
                    st[ns].Score=score;
                    st[ns].sclass= klass;
                    st[ns].type= type;
                    st[ns].subtype= subtype;
                    ns++;
                }

                is++;
                s= s.substring(1);  // s++
            }
        }
        /* case C: Using Markov chain with order > 1 */
        else
        {
            /* discovering splice sites with current profile */
            while (s.length()< (p.dimension- 1) && (is < right- left + 1) && (ns<nsites))   // *(s+p.dimension-1) &&
            {
                /* is = 0..right */
                score=0f;
                for (int i=0;i<p.dimension;i++)
                {
                    /* i is the position inside the region */
                    /* 5 is used because there are A,C,G,T and N */
                    index = OligoToInt(s.substring(i- p.order), p.order+1,5);  // s + i - p.order

                    if (index >= p.dimensionTrans)
                        score = score + -GeneIDconstants.INFI;
                    else
                        score = score + p.transitionValues[i][index];
                }
                score = p.afactor + (p.bfactor * score);
                if(GeneIDsettings.UTR!= 0){
                    score = score - peakEdgeScore(left + is + p.offset, strand, external, l1, l2, 6);
                }
                if (score >= p.cutoff)
                {
                    st[ns].Position = left + is + p.offset;
                    st[ns].Score=score;
                    st[ns].sclass= klass;
                    st[ns].type= type;
                    st[ns].subtype= subtype;
                    ns++;
                }
                is++;
                s= s.substring(1);  // s++
            }
        }

        if (ns >= nsites)
            Log.error("Too many predicted sites: decrease RSITES parameter");

        return(ns);
    }


    float peakEdgeScore(int Position,
                        int Strand,
                        PackExternalInformation external,
                        int l1, long l2, int win)
    {
        int index;
        short trueFrame;
        int relPos;
        float Score = 0;
        float factor = 10;
/*   int win = 6; */
/*   char mess[MAXSTRING]; */

        relPos = Position - l1 + GeneIDconstants.COFFSET;

        if (Strand == GeneIDconstants.FORWARD)
            index = 0;
        else
            index = GeneIDconstants.FRAMES;

        /* Frame blast definition: according to the sequence start */
        trueFrame = 0;
        trueFrame += index;

        /* Access the sr array to obtain the homology score for current score */
        if (((relPos -win) >=0) && ((relPos +win) < LENGTHSi) && ((Position+win)<l2)){
            Score = (factor*(
                    ((external.sr[trueFrame][relPos + win] - external.sr[trueFrame][relPos])/win)
                            - ((external.sr[trueFrame][relPos] - external.sr[trueFrame][relPos -win])/win)
            )
            );
/*     sprintf(mess,"relPos: %ld\nslope 2: %f\nslope 1: %f",relPos,((external.sr[trueFrame][relPos +win] - external.sr[trueFrame][relPos])/win),((external.sr[trueFrame][relPos] - external.sr[trueFrame][relPos -win])/win)); */
/*     printMess(mess); */
/*     sprintf(mess,"peakEdgeScore: %f",(factor*(((external.sr[trueFrame][relPos +win] - external.sr[trueFrame][relPos])/win)  */
/* 					    - ((external.sr[trueFrame][relPos] - external.sr[trueFrame][relPos -win])/win)))); */
/*     printMess(mess); */

        }

        return(Score);
    }

    double computeU2BranchProfile(String s, int positionAcc, long limitRight, Profile p, Site splicesite) {

        double maxScore = -GeneIDconstants.INF;

        /*      char mess[MAXSTRING];  */
        int end = (int) Math.min(positionAcc - p.dist + p.dimension - p.offset,limitRight);
        end = (int) Math.min(end, positionAcc);
        int Opt = (int) Math.max(0,positionAcc - p.opt_dist);
        for (int i= (int) Math.max(p.order,positionAcc - p.acc_context);
             i + p.dimension <= end;
             i++) {
            /* Applying the additional profile */
            float score= 0f;
            for (int j= 0; j < p.dimension; j++) {
                /* i is the position inside the region */
                /* 5 is used because there are A,C,G,T and N */
                int index = OligoToInt(s.substring(i + j - p.order), p.order+1,5);

                if (index >= p.dimensionTrans)
                    score = score + -GeneIDconstants.INFI;
                else
                    score = score + p.transitionValues[j][index];
            }

            score = score - p.penalty_factor * (((float) (Math.abs(i + p.offset -Opt))/((float)(p.acc_context - p.offset - p.opt_dist))) *
                    ((float)(Math.abs(i + p.offset- Opt))/((float)(p.acc_context - p.offset - p.opt_dist))));
            score = p.afactor + (p.bfactor * score);
            if ((score >= maxScore)&&(score > p.cutoff)){
                maxScore = score;
                splicesite.PositionBP = i + p.offset - positionAcc;
            }
        }

        /* Cutoff for BranchPoint and PPtracts are useless */
/*   if (maxScore < p.cutoff) */
/*   	maxScore = 0f; */

        return maxScore;
    }


    double computeU2PPTProfile(String s,
                              int positionAcc,
                              long limitRight,
                              Profile p,
                              Site splicesite) {

        double maxScore = -GeneIDconstants.INF;

        int end = (int) Math.min(positionAcc,limitRight);
        for (int i = (int) Math.max(p.order, positionAcc - p.dist - p.dimension);
             i + p.dimension <= end;
             i++) {
            /* Applying the additional profile */
            float score= 0f;
            for (int j= 0; j < p.dimension; j++) {
                /* i is the position inside the region */
                /* 5 is used because there are A,C,G,T and N */
                int index = OligoToInt(s.substring(i + j - p.order), p.order+ 1,5);

                if (index >= p.dimensionTrans)
                    score = score + -GeneIDconstants.INFI;
                else
                    score = score + p.transitionValues[j][index];
            }

            if ((score >= maxScore)&&(score > p.cutoff)){
                maxScore = score;
                splicesite.PositionPPT = i + p.offset - positionAcc;
            }
        }

        /* Cutoff for BranchPoint and PPtracts are useless */
/*   if (maxScore < p.cutoff) */
/*   	maxScore = 0f; */

        return maxScore;
    }


}
