package CoInfecProfiler;

import java.time.format.DateTimeFormatter;
import java.util.HashSet;

public class Main {
	
	public static void main(final String[] args) throws Exception {
		System.out.println("Welcome to use CoInfecProfiler! If you have any questions, please refer to "
		 		+ "caochen@njmu.edu.cn");
		final String[] supported_functions_array = 
			{"SARS-CoV-2-Reconstruction-VCF"};
		final HashSet<String> supported_functions = new HashSet<String>();
        for (int k = 0; k < supported_functions_array.length; ++k) {
            supported_functions.add(supported_functions_array[k]);
        }
        final String function = args[0];
        
        if (!supported_functions.contains(function)) {
            System.out.println("Function " + function + " is not supported. A typo?");
            System.exit(0); 
        }
        
		if (function.equals("SARS-CoV-2-Reconstruction-VCF")) {
        	final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
        	String config_path = args[1];
        	final SARS_CoV_2 sc = new
        			SARS_CoV_2(config_path);
			sc.run_bwa();
			sc.run_GATK();
        	sc.readvcf();
        	sc.variantmaker(); 
        	sc.vefgenerater();
        	sc.readref(); 
        	sc.readmatrix();
        	sc.l0l1init();
        	sc.rscript();
        	sc.run_R();
        	sc.userformat();
        	System.exit(0);
        }
	}
}

