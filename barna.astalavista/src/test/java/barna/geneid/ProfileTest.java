package barna.geneid;

import barna.commons.Execute;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Created with IntelliJ IDEA.
 * User: micha
 * Date: 9/8/12
 * Time: 6:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class ProfileTest {

    @BeforeClass
    public static void initExecuter() {
        Execute.initialize(2);
    }

    @AfterClass
    public static void shutdownExecuter() {
        Execute.shutdown();
    }

    @Test
    public void testReadProfile() throws Exception {
        GeneIDsettings settings= new GeneIDsettings();
        GParam[] isochores= Profile.readParam(GeneIDconstants.PARAMETERFILE, settings);
    }

    @Test
    public void testBuildDonors() throws Exception {
        GeneIDsettings settings= new GeneIDsettings();
        GParam[] isochores= Profile.readParam(GeneIDconstants.PARAMETERFILE, settings);
        System.currentTimeMillis();

        GeneID myGID= new GeneID();
        myGID.scoreDonor("GTACCCC", isochores[0].DonorProfile);
//        myGID.buildDonors(
//                "",
//                (short) 0,
//                "", // type
//                "", // sub-type
//                isochores[0].DonorProfile, // Profile
//                new Site[1], //Site[] st
//                0,  // l1
//                0,  // l2
//                0,  // ns
//                1,  // nsites
//                GeneIDconstants.FORWARD,    // Strand
//                null    // Pack external
//        );
    }

}
