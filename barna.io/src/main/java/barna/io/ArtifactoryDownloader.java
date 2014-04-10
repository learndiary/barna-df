package barna.io;

import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * Class to download stuff from the artifactory.
 */
public class ArtifactoryDownloader {

    /**
     * The test data target directory
     */
    private static String targetDirectory;

    /**
     * Complete Gencode annotation
     */
    private static File gencodeFile= getGencodeFile();


    public static File getGencodeFile() {
        if (gencodeFile== null)
            gencodeFile= prepareTestData();
        return gencodeFile;
    }

    static {
        // initialize defaults
        //artifactoryUrl = System.getProperty("testdata.artifactory", "http://sammeth.net/artifactory")
        //repository = System.getProperty("testdata.repository", "repo")
        //artifact = System.getProperty("testdata.artifact", "barna/test_data-1.0.zip")
        targetDirectory = System.getProperty("testdata.target", new File("").getAbsolutePath());
    }

    private static File prepareTestData() {

        System.out.println("Checking for test data");
        //JsonSlurper slurper = new JsonSlurper();
        URL dataFQN = null;
        try {
            dataFQN = new URL(
                    // ${artifactoryUrl}/api/storage/${repository}/${artifact}
                    " http://sammeth.net/artifactory/repo/gencode_v12_gtf/gencode_v12_gtf/12/gencode_v12_gtf-12.gz");
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }
        String[] tokens= dataFQN.getFile().split("/");
        String fileName = tokens[tokens.length- 1];


        File targetFile = new File(targetDirectory, fileName);
        File test_data_dir = new File(targetDirectory, "test_data");
        File md5File = new File(targetDirectory, "${fileName}.md5");

        //metaData = slurper.parseText(dataFQN.openStream().text);
        //def md5sum = metaData.checksums['md5']

        boolean invalid_file = false;
        if (!targetFile.exists()) {
            // || !md5File.exists()) || md5File.readLines()[0].trim() != md5sum) {
            invalid_file = true;
        }

        // download
        if (invalid_file) {
            System.out.println("Downloading test data from ${dataFQN.toExternalForm()}");
            if (targetFile.exists()) {
                targetFile.delete();
            }
            if (test_data_dir.exists()) {
                FileHelper.rmDir(new File(targetDirectory, "test_data"));
            }
            if (md5File.exists()) {
                md5File.delete();
            }
            // download
            //out << new URL(metaData['downloadUri']).openStream();
            try {
                InputStream in= dataFQN.openStream();
                OutputStream out = new FileOutputStream(targetFile);
                byte[] buffer = new byte[1024];
                int len;
                while ((len = in.read(buffer)) != -1) {
                    out.write(buffer, 0, len);
                }
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            // compare md5
            //if (md5sum != generateMD5(targetFile)) {
            //    throw new RuntimeException("Test Data downlaoded but md5 sums do not match !");
            //}

            // write md5 file
            //md5File.write(md5sum)
        }

        File gencodeFile= new File(FileHelper.stripExtension(targetFile.getAbsolutePath()));
        if (!gencodeFile.exists()) {
            System.out.println("Unzipping test data");
            //FluxCapacitorRunner.unzip(targetFile, targetDirectory);
            try {
                FileHelper.inflate(targetFile, gencodeFile, FileHelper.getCompression(targetFile));
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Test data available");
        return gencodeFile;
    }

}
