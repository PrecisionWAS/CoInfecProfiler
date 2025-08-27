package CoInfecProfiler;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;


public class VariantsCall {
	public String gatk;
	public String bwa;
	public String samtools;
	public String project_name;
	public String fq_1;
	public String fq_2;
	public String reference;
	public String out_dir;
	public String java;
	public int min_read_cov= 20;
	public String python;
	public String spades;
	public String ref;
	public  HashMap<Integer, String> Pos2Ref;
	public  HashMap<Integer, String> Pos2Alt;
	public  ArrayList<Integer> PosArr;
	
	
	public VariantsCall(String parameter_file) throws IOException {
		InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
//        this.gatk = prop.getProperty("GATK") ;
        this.bwa = prop.getProperty("bwa") ;
        this.samtools = prop.getProperty("samtools") ;
        this.fq_1 = prop.getProperty("Fastq_1_Path") ;
        this.fq_2 = prop.getProperty("Fastq_2_Path") ;
        this.out_dir= prop.getProperty("Output_Path")+"/" ;
        this.reference= prop.getProperty("Reference_Seq") ;
        this.project_name = prop.getProperty("Proj_Name") ;
        is.close();
	}
	
	public VariantsCall(String parameter_file, int meta) throws IOException {
		InputStream is = new FileInputStream(parameter_file);
        Properties prop = new Properties();
        prop.load(is);
        this.bwa = prop.getProperty("bwa") ;
        this.samtools = prop.getProperty("samtools") ;
        this.fq_1 = prop.getProperty("Fastq_1_Path") ;
        this.fq_2 = prop.getProperty("Fastq_2_Path") ;
        this.out_dir= prop.getProperty("Output_Path")+"/" ;
        this.reference= prop.getProperty("Reference_Seq") ;
        this.project_name = prop.getProperty("Proj_Name") ;
        this.python = prop.getProperty("PYTHON") ;
        this.spades = prop.getProperty("Spades") ;
        is.close();
	}
	
	public void process_rawhaplo()  throws IOException, InterruptedException {
		
		BufferedReader br = new BufferedReader(new FileReader(String.valueOf(this.out_dir)+ 
				"/intermediate/metaSPAdes/contigs.raw.haps"));

		ArrayList<Integer >  pos_set = new ArrayList<Integer>();
		
		String line= "";
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
			for (int i=0;i< tmp.length ;i++) {
				if (i>1) {
					pos_set.add(Integer.parseInt(tmp[i]) );
				}
			}
		}
		br.close();
		
		BufferedReader br2 = new BufferedReader(new FileReader(String.valueOf(this.out_dir)+ 
				"/intermediate/metaSPAdes/contigs.raw.haps"));
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				String.valueOf(this.out_dir)+ "/intermediate/metaSPAdes/contigs.haps", false));
				
		while ((line = br2.readLine()) != null) {
			ArrayList<Integer >  pos_arr = new ArrayList<Integer>();
			ArrayList<String >   flag_arr = new ArrayList<String>();
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
			
			int start =  Integer.parseInt(tmp[0]);
			int end =  Integer.parseInt(tmp[1]);
			
			for (int i=0;i< tmp.length ;i++) {
				if (i>1) {
					pos_arr.add(Integer.parseInt(tmp[i]));
					flag_arr.add(tmp[i]+"=1"  );
				}	
			}
			
			for (int i=0;i< pos_set.size() ;i++) {
				if ( ( pos_set.get(i) > start )&&  
						( pos_set.get(i) < end ))  {
					boolean flag= true;
					for (int j=0;j< pos_arr.size(); j++) {
//						System.out.println(pos_arr.get(j)+"\t"+ pos_set.get(i));
//						System.out.println(pos_set.get(i)-pos_arr.get(j));
						if( ( pos_set.get(i)-pos_arr.get(j) )==0 )    {
//							System.out.println(pos_arr.get(j)+"\t"+ pos_set.get(i));
							flag=false;
						}
					}
					if (flag==true) { 
						pos_arr.add(pos_set.get(i) );
						flag_arr.add(Integer.toString(pos_set.get(i) ) +"=0"  );
					}
				}	
			}
			
			for (int i=0; i< pos_arr.size(); i++) {
				for (int j= (i+1) ; j < pos_arr.size(); j++) {
					if (pos_arr.get(i) >  pos_arr.get(j)) {
					
						int tmp_int = pos_arr.get(i);
						pos_arr.set(i , pos_arr.get(j));
						pos_arr.set(j , tmp_int);
						
						String tmp_str = flag_arr.get(i);
						flag_arr.set(i, flag_arr.get(j));
						flag_arr.set(j,  tmp_str);
					}
				}
			}
							
			String ss= tmp[0]+"\t"+ tmp[1];
			for (int i=0; i< flag_arr.size();i ++) {
				ss =ss +"\t" + flag_arr.get(i);
			}
			bw.write(ss+"\n");
		}
		bw.close();
		br2.close();	
	}
	
	
	
	public void read_hap()  throws IOException, InterruptedException {
		BufferedReader br = new BufferedReader(new FileReader(String.valueOf(this.out_dir)+ 
				"/intermediate/metaSPAdes/contigs.raw.haps"));
		ArrayList<Integer >  pos_set = new ArrayList<Integer>();
		String line= "";
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
			for (int i=0;i< tmp.length ;i++) {
				if (i>1) {
					boolean flag =true;
					for (int j=0; j< pos_set.size();j++) {
						if ( pos_set.get(j)== Integer.parseInt(tmp[i]) ) {
							flag= false;
						}
					}
					if (flag) {
						pos_set.add(Integer.parseInt(tmp[i]) );
					}
				}
			}
		}
		br.close();
		for (int i=0;i< pos_set.size();i++) {
			for (int j=(i+1);j< pos_set.size();j++) {
				if  ( pos_set.get(i) > pos_set.get(j) ) {
					int tmp_int = pos_set.get(i);
					pos_set.set(i , pos_set.get(j));
					pos_set.set(j , tmp_int);
				}
			}
		}
		
		HashMap<Integer, String> Pos2SNP = new HashMap<Integer, String>();
		
		this.PosArr = new ArrayList<Integer>();
		for (int i=0;i< pos_set.size();i++ ) {
			this.PosArr.add(pos_set.get(i)); 
			Pos2SNP.put(pos_set.get(i) , "");
		}
		
		BufferedReader br2 = new BufferedReader(new FileReader(String.valueOf(this.out_dir)+ 
				"/intermediate/metaSPAdes/contigs.haps"));	
		int line_index =0;
		while ((line = br2.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
			for (int i=0;i< tmp.length ;i++) {
				if (i>1) {
					String [] tmp2= tmp[i].split("=");
					int pos =  Integer.parseInt(tmp2[0]);
					Pos2SNP.put(pos, Pos2SNP.get(pos)+ 
							Integer.toString(line_index ) +":"+ tmp2[1]+";");
				}	
			}
			line_index++;	
		}
		br2.close();
		
		ArrayList<String >  Ref_Arr = new ArrayList<String>();
		ArrayList<String >  Pos_Arr = new ArrayList<String>();
		ArrayList<String >  Hap = new ArrayList<String>();
		
		Ref_Arr.add("");
		Pos_Arr.add("");
		Hap.add("-1");
		for (int i=0;i< pos_set.size(); i++) {
			String str = Pos2SNP.get(pos_set.get(i ));
			String tmp[] = str.split(";");
			if ( tmp.length>1 ) {
				ArrayList<String >  Tmp_Ref_Arr = new ArrayList<String>();
				ArrayList<String >  Tmp_Pos_Arr = new ArrayList<String>();
				ArrayList<String >  Tmp_Hap = new ArrayList<String>();
				
				for (int j=0;j < Ref_Arr.size();j++) {
					Tmp_Ref_Arr.add(Ref_Arr.get(j));
					Tmp_Pos_Arr.add(Pos_Arr.get(j));
				}
				
				for (int j=0;j< Hap.size();j++) {
					Tmp_Hap.add(Hap.get(j));
				}
				
				Ref_Arr.clear();
				Pos_Arr.clear();
				Hap.clear();
				boolean flag = false;
				for (int j=0;j< tmp.length  ;j++) {
					String h= tmp[j].split(":")[0];
					String allele = tmp[j].split(":")[1];
					for (int k=0; k< Tmp_Hap.size();k ++) {
						if (Tmp_Hap.get(k).equals( h)) {
							flag= true;
						}
					}
				}
				
				if (flag==true) {
					for (int j=0; j< Tmp_Hap.size(); j++) {
						for (int k=0;k< tmp.length  ;k++) {
							String h= tmp[k].split(":")[0];
							String allele = tmp[k].split(":")[1];
							if (Tmp_Hap.get(j).equals( h)  ) {
								Ref_Arr.add( Tmp_Ref_Arr.get(j)+ allele );
								Pos_Arr.add(Tmp_Pos_Arr.get(j)+  Integer.toString(pos_set.get(i )) +";" ) ;
								Hap.add(h );
							}
						}
					}
				} else {
					for (int j=0; j< Tmp_Hap.size(); j++) {
						for (int k=0;k< tmp.length  ;k++) {
							String h= tmp[k].split(":")[0];
							String allele = tmp[k].split(":")[1];
							if (!Tmp_Hap.get(j).equals( h)  ) {
								Ref_Arr.add( Tmp_Ref_Arr.get(j)+ allele );
								Pos_Arr.add(Tmp_Pos_Arr.get(j)+  Integer.toString(pos_set.get(i )) +";" ) ;
								Hap.add(h );
							}
						}
					}
				}
			}
		}
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				String.valueOf(this.out_dir)+ 
				"/intermediate/metaSPAdes/contigs.SNV", false));
		
		if (Pos_Arr.size() > 0) {
			System.out.println(Pos_Arr.get(0));
			String tmp[] = Pos_Arr.get(0).split(";");
			for (int i=0; i< tmp.length ; i++) {
				System.out.println(tmp[i] );
				int pos= Integer.parseInt ( tmp[i]) ;
				String ref = this.Pos2Ref.get(pos);
				String alt= this.Pos2Alt.get(pos );
				String ss= Integer.toString( pos)+"\t"+ ref+"\t"+alt;
				bw.write(ss+"\n");
			}
		}
		
		
		bw.close();
		
		new File(this.out_dir+"/intermediate/1_"+ Integer.toString(this.ref.length() )).mkdir();
		
		BufferedWriter bw2 = new BufferedWriter(new FileWriter(
				this.out_dir+"/intermediate/1_"+ Integer.toString(this.ref.length() )+ "/"+this.project_name+
				"_inter.freq", false));
		String ss="Hap_ID";
		for (int i=0;i< Ref_Arr.size();i++) {
			ss=ss+"\th"+Integer.toString(i);
		}
		bw2.write(ss+"\n");
		ss="Freq";
		for (int i=0;i< Ref_Arr.size();i++) {
			double n = Ref_Arr.size();
			ss=ss+"\t" + Double.toString( 1/ n);
		}
		bw2.write(ss+"\n");
		
		if ( Pos_Arr.size()> 0 ) {
			for (int i=0; i< Ref_Arr.get(0).length();i++) {
				String tmp[] = Pos_Arr.get(0).split(";"); 
				ss="0;"+tmp[i]+";"+tmp[i]+";0:1";
				
				for (int j=0;j< Ref_Arr.size();j++) {
					ss=ss+"\t"+ Ref_Arr.get(j).substring(i, i+1);
				}
				bw2.write(ss+"\n");
			}
		}
		
		bw2.close();
		
//		for (int i=0;i< Ref_Arr.size(); i++) {
//			System.out.println(Ref_Arr.get(i)+"\t"+ Pos_Arr.get(i)+"\t"+Hap.get(i)  );
//		}
	}
		
	public void sam_generate_haplo()  throws IOException, InterruptedException {
		this.Pos2Ref = new HashMap<Integer, String>  ();
		this.Pos2Alt = new HashMap<Integer, String>  ();
		BufferedReader br = new BufferedReader(new FileReader(this.reference));
		String line= "";
		this.ref="";
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			if ( !line.substring(0, 1).equals(">")) {
				this.ref =this.ref +line.toUpperCase();
			}	
		}
		br.close();
		String sam_fil = String.valueOf(this.out_dir)+ 
				"/intermediate/metaSPAdes/contigs.sam";
		BufferedReader br_sam = new BufferedReader(new FileReader(sam_fil));
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				String.valueOf(this.out_dir)+ "/intermediate/metaSPAdes/contigs.raw.haps", false));
		
		final HashSet<String> cigar_symbols = new HashSet<String>();
		cigar_symbols.add("M");		
		cigar_symbols.add("D");		
		cigar_symbols.add("I");		
		cigar_symbols.add("H");
		cigar_symbols.add("S");
		while ((line = br_sam.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp = line.split("\t");
//			System.out.println(line);
			if ( ( !line.substring(0, 1).equals("@")) && ( !tmp[3].equals("*")) ){
				String ss="";
				int Start_Pos= Integer.parseInt(tmp[3]);
				String read= tmp[9];
				String cigar= tmp[5];
//				System.out.println(cigar);
				int cigar_ps =0;
				int read_ps=0;
				
				int ref_ps= Start_Pos;
				int p=0;
				while (p < (cigar.length()-1)  ) {
					p++;
					if ( cigar_symbols.contains( cigar.substring(p, p+1) ) ){
						int len= Integer.parseInt(cigar.substring( cigar_ps, p));
//						System.out.println(Integer.toString(len)+ " "+ cigar.substring(p, p+1));
						cigar_ps=p+1;
						if ( cigar.substring(p, p+1).equals("M") ) {
//							System.out.println(  read.substring(read_ps, read_ps+len));
//							System.out.println(  this.ref.substring(ref_ps-1, ref_ps-1+len));
							for (int i=0; i< len; i++) {
								if ( !read.substring(read_ps, read_ps+1).equals( 
										this.ref.substring(ref_ps-1, ref_ps))) {
									ss= ss+ "\t"+ Integer.toString(ref_ps);
									this.Pos2Ref.put(ref_ps , this.ref.substring(ref_ps-1, ref_ps)); 
									this.Pos2Alt.put(ref_ps , read.substring(read_ps, read_ps+1)); 
								}
								read_ps++;
								ref_ps++;
							}
						} else if ( cigar.substring(p, p+1).equals("I") ) {
							for (int i=0; i< len; i++) {
								read_ps=read_ps+1;
							}
						} else if ( cigar.substring(p, p+1).equals("D") ) {
							for (int i=0; i< len; i++) {
								ref_ps=ref_ps+1;
							}
						}	else if ( cigar.substring(p, p+1).equals("S") ) {
							for (int i=0; i< len; i++) {
								read_ps=read_ps+1;
							}
						}
					} 
				}		
				bw.write(tmp[3]+"\t"+ Integer.toString( ref_ps) + ss + "\n"); 
			}	
		}
		bw.close();
		br_sam.close();
	}

	
	public void refindex()  throws IOException, InterruptedException {
		new File(String.valueOf(this.out_dir) ).mkdir();
		//$bwa mem $ref $prefix\.bwa.read1.fastq $prefix\.bwa.read2.fastq 
//		Runtime r=Runtime.getRuntime();
//		String cmd = this.bwa +" mem "+this.reference+" "+this.fq_1+" "+ this.fq_2+" > "+this.out_dir+"tmp.sam";
//		System.out.println(cmd);
//		Process p=r.exec(cmd);
//		InputStream is=p.getInputStream();
//		InputStreamReader ir=new InputStreamReader(is);
//		BufferedReader br=new BufferedReader(ir);
//		String str=null;
//		while((str=br.readLine())!=null){
//			System.out.println(str);
//		}
//	    int ret=p.waitFor();
//	    int exit_v=p.exitValue();
//	    System.out.println("return value:"+ret);
//	    System.out.println("exit value:"+exit_v);
		new File(String.valueOf(this.out_dir) ).mkdir();
		new File(String.valueOf(this.out_dir)+"/intermediate" ).mkdir();
	    ProcessBuilder CMDLine = new ProcessBuilder(this.bwa,
	                "index", 
	                this.reference
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
	            String line;
	            while ((line = br.readLine()) != null) {
	            	line= line.replace("\n", "").replace("\r", "");
	               
	            }
	            CMDProcess.waitFor();
	            System.out.println("Finished indexing the reference.");
	            br.close();
	}
	
	public void metaSPAdes()  throws IOException, InterruptedException {
		new File(String.valueOf(this.out_dir) ).mkdir();
		new File(String.valueOf(this.out_dir)+"/intermediate" ).mkdir();
		new File(String.valueOf(this.out_dir)+"/intermediate/metaSPAdes" ).mkdir();
		
	    ProcessBuilder CMDLine = new ProcessBuilder(this.python,
	    		this.spades,
                "-o",
                String.valueOf(this.out_dir)+"/intermediate/metaSPAdes",
                "--meta",
                "-1",
                this.fq_1,
                "-2",
                this.fq_2
               );
	            Process CMDProcess = CMDLine.start();
	            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
	            String line;
	            while ((line = br.readLine()) != null) {
	            	line= line.replace("\n", "").replace("\r", "");
	            }
	            CMDProcess.waitFor();
	            System.out.println("Finished running metaSPAdes.");
	            br.close();   
	}
	
	public void bwa2()  throws IOException, InterruptedException {
		new File(String.valueOf(this.out_dir) ).mkdir();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(String.valueOf(this.out_dir)+
				"/intermediate/metaSPAdes/contigs.sam", false));
		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.bwa,
	                "mem", 
	                this.reference, 
	                String.valueOf(this.out_dir)+"/intermediate/metaSPAdes/contigs.fastq"
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
	            String line;
	            while ((line = br.readLine()) != null) {
	            	line= line.replace("\n", "").replace("\r", "");
	                bw.write(line+"\n");
	            }
	            CMDProcess.waitFor();
	            System.out.println("Finished aligning metaSPAdes reads.");
	            br.close();
		bw.close();
	}
	
	
	public void contig2fastq()  throws IOException, InterruptedException {
		
		ArrayList<String >  reads = new ArrayList<String>();
		ArrayList<String >  covs = new ArrayList<String>();
		
		BufferedReader br = new BufferedReader(new FileReader(String.valueOf(this.out_dir)+
				"/intermediate/metaSPAdes/contigs.fasta"));
		
		BufferedWriter bw = new BufferedWriter
				(new FileWriter(String.valueOf(this.out_dir)+"/intermediate/metaSPAdes/contigs.fastq", 
						false));
		
		String line ="";
		String contig ="";
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			if  ( line.substring(0, 1).equals(">") ) {
				String[] tmp = line.split("_cov_");
				covs.add(tmp[1]);
//				System.out.println(contig.length());
				if ( contig.length()> 500 ) {
					reads.add(contig);
//					System.out.println(contig);
				}
				contig ="";
			}else {
				contig = contig + line;
			}
		}
		br.close();
		
		for (int i=0;i< reads.size();i++ ) {
			String ss="@CONTIG_"+ Integer.toString(i)+"_"+ covs.get(i);
			bw.write(ss+"\n");
			bw.write(reads.get(i)+"\n");
			ss="+_"+ Integer.toString(i)+"_"+ covs.get(i);
			bw.write(ss+"\n");
			ss="";
			for (int j=0; j< reads.get(i).length();j++) {
				ss=ss+'H';
			}
			bw.write(ss+"\n");
		}
		bw.close();
	}
	
	public void bwa()  throws IOException, InterruptedException {
		new File(String.valueOf(this.out_dir) ).mkdir();
		//$bwa mem $ref $prefix\.bwa.read1.fastq $prefix\.bwa.read2.fastq 
//		Runtime r=Runtime.getRuntime();
//		String cmd = this.bwa +" mem "+this.reference+" "+this.fq_1+" "+ this.fq_2+" > "+this.out_dir+"tmp.sam";
//		System.out.println(cmd);
//		Process p=r.exec(cmd);
//		InputStream is=p.getInputStream();
//		InputStreamReader ir=new InputStreamReader(is);
//		BufferedReader br=new BufferedReader(ir);
//		String str=null;
//		while((str=br.readLine())!=null){
//			System.out.println(str);
//		}
//	    int ret=p.waitFor();
//	    int exit_v=p.exitValue();
//	    System.out.println("return value:"+ret);
//	    System.out.println("exit value:"+exit_v);
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(this.out_dir+"/intermediate/"+ 
				this.project_name+".sam", false));
		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.bwa,
	                "mem", 
	                this.reference, 
	                this.fq_1,
	                this.fq_2
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
	            String line;
	            while ((line = br.readLine()) != null) {
	            	line= line.replace("\n", "").replace("\r", "");
	                bw.write(line+"\n");
	            }
	            CMDProcess.waitFor();
	            System.out.println("Finished aligning reads.");
	            br.close();
		bw.close();
	}
	
	public void sam2bam()  throws IOException, InterruptedException {

		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.samtools,
	                "view", 
	                "-Shub", 
	                this.out_dir+"/intermediate/"+this.project_name+".sam",
	                "-o",
	                this.out_dir+"/intermediate/"+this.project_name+".raw.bam"
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            
	            CMDProcess.waitFor();
	            System.out.println("Finished Converting sam to bam.");
	}
	
	public void sort()  throws IOException, InterruptedException {

		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.samtools,
	                "sort", 
	                "-o", 
	                this.out_dir+"/intermediate/"+this.project_name+".bam",
	                this.out_dir+"/intermediate/"+this.project_name+".raw.bam"
	                
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            
	            CMDProcess.waitFor();
	            System.out.println("Finished sorting the bam file.");
	            
	   
	}
	
	public void index()  throws IOException, InterruptedException {

		
	    
	    ProcessBuilder CMDLine = new ProcessBuilder(this.samtools,
	                "index", 
	                this.out_dir+"/intermediate/"+this.project_name+".bam"
	               );
	            Process CMDProcess = CMDLine.start();
	            
	            CMDProcess.waitFor();
	            System.out.println("Finished index the bam file.");
	}
	
	
	
	public void AddOrReplaceReadGroups()  throws IOException, InterruptedException {

	    ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "AddOrReplaceReadGroups",
	                "-I",
	                this.out_dir+this.project_name+".srt.bam",
	                "-O",
	                this.out_dir+this.project_name+".rg.bam",
	                "-R",
	                this.reference,
	                "-ID",
	                this.project_name,
	                //-LB NPD -PL Illumina -PU NPD -SM 
	                "-LB",
	                "NPD",
	                "-PL",
	                "Illumina",
	                "-PU",
	                "NPD",
	                "-SM",
	                this.project_name
	                
	                
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
	            System.out.println("Finished AddOrReplaceReadGroups the bam file.");
	}	
	
	
//	public void call()  throws IOException, InterruptedException {
//
//	    ProcessBuilder CMDLine = new ProcessBuilder(this.java,
//	                "-jar", 
//	                this.gatk,
//	                "HaplotypeCaller",
//	                "-I",
//	                this.out_dir+this.project_name+".rg.bam",
//	                "-R",
//	                this.reference,
//	             
//	                "-O",
//	                this.out_dir+this.project_name+".tmp.vcf"
//
//	               );
//	            Process CMDProcess = CMDLine.start();
//	            CMDProcess.waitFor();
//	            System.out.println("Finished calling the variants using GATK.");
//	}	
	
	public void call()  throws IOException, InterruptedException {
//		$java -jar $gatk  HaplotypeCaller -R  $ref -I $inbam -ERC  GVCF -ploidy 8 --heterozygosity 0.01  --max-alternate-alleles 1 -O $outgvcf
	    ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "HaplotypeCaller",
	                "-I",
	                this.out_dir+this.project_name+".rg.bam",
	                "-R",
	                this.reference,
	                "-ERC",
	                "GVCF",
	                "-ploidy",
	                "8",
	                "--heterozygosity",
	                "0.01",
	                "-max-alternate-alleles",
	                "1",
	                "-O",
	                this.out_dir+this.project_name+".raw.g.vcf"
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
	            System.out.println("Finished calling the variants using GATK.");
	}	
	
	
	public void GenotypeGVCFs()  throws IOException, InterruptedException {
//		$java -jar $gatk   GenotypeGVCFs -R $ref -V  $prefix\.g.vcf -ploidy 8 -O $prefix\.raw.vcf	    
		ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "GenotypeGVCFs",
	                "-R",
	                this.reference,
	                "-ploidy",
	                "8",
	                "-V",
	                this.out_dir+this.project_name+".raw.g.vcf",
	                "-O",
	                this.out_dir+this.project_name+".raw.vcf"
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
//	            System.out.println("Finished calling the variants using GATK.");
	}	
	
	public void SelectVariants()  throws IOException, InterruptedException {
//		$java -jar $gatk  SelectVariants -R $ref -V  $prefix\.raw.vcf -O $prefix_vcf\.vcf    
		ProcessBuilder CMDLine = new ProcessBuilder(this.java,
	                "-jar", 
	                this.gatk,
	                "SelectVariants",
	                "-R",
	                this.reference,

	                "-V",
	                this.out_dir+this.project_name+".raw.vcf",
	                "-O",
	                this.out_dir+this.project_name+".tmp.vcf"
	               );
	            Process CMDProcess = CMDLine.start();
	            CMDProcess.waitFor();
//	            System.out.println("Finished calling the variants using GATK.");
	}	
	
	
	public double  freq_get(String ss)  throws IOException, InterruptedException {
		
		String[] tmp = ss.split("\t");  
		String tt = tmp[9];
		String[] tmp_tt = tt.split(":");
		String yy= tmp_tt[1];
		double ref =Double.parseDouble(yy.split(",")[0]);
		double alt =Double.parseDouble(yy.split(",")[1]);
		if ((ref+alt)<min_read_cov) {
			return 0.0;
		} else {
			return alt/(ref+alt);
		}
			
		
	}
	public void filter_low()  throws IOException, InterruptedException {
		int snp_end= -1;
		BufferedWriter bw = new BufferedWriter
				(new FileWriter(this.out_dir+this.project_name+".vcf", false));
		BufferedReader br = new BufferedReader(new FileReader(this.out_dir+this.project_name+".tmp.vcf" ));
		String line;
		while ((line = br.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			if (line.substring(0, 1).equals("#")) {
				bw.write(line+"\n");
			}else {
				double freq = freq_get(line);
				String []tmp =line.split("\t"); 
				if ((freq>0.01) &&  (freq<0.99) && (tmp[3].length()==1) && 
						(tmp[4].length()==1)) {
					int this_start = Integer.parseInt(tmp[1]);
					int this_end = Integer.parseInt(tmp[1]) + tmp[3].length()-1;
					if (this_start> snp_end) {
						bw.write(line+"\n");
						snp_end= this_end;
					}
				}
			}
		}
		br.close();
		bw.close();
	}	

}
