import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;

/*
 * Here we use the dataset from the internet,which is not suited to data processing,
 * so we use this code to convert the format of the dataset.
 */
public class preTreatment {

	
 @SuppressWarnings("unchecked")
 public static void main(String[] args) {
	
	
  int c = Integer.valueOf(args[0]);		//columns number of input file
  String f = args[1];					//the path of input file
  String fout = args[2];				//the path of output file
	
	
	
  try {
	  long startTime = System.currentTimeMillis();
	  File file = new File(fout);

	  //If file doesn't exists,create it.
	  if (!file.exists()) {
		  file.createNewFile();
	  }

	  try {
		  BufferedReader reader = new BufferedReader(new FileReader(f));
          String line = null;
          @SuppressWarnings("rawtypes")
          Vector Data[] = new Vector[c];
          System.out.println("Start.");
          FileWriter fw = new FileWriter(file.getAbsoluteFile());
    	  BufferedWriter bw = new BufferedWriter(fw);
    	  int j = 0;
          while((line=reader.readLine())!=null){
        	    
              	String item[] = line.split("  ");
				//find relationships  between genes.  
				for(int i = 0;i < item.length;i++){
					
					if(Data[j] == null){
              			Data[j] = new Vector<String>();
              		}
					
					Data[j].addAll(Arrays.asList( "V" +  (i+1) +"*" + item[i]));
				}
				//System.out.println(Data[j]);
				j++;
          	}
          for(int i = 0;i < c ;i++){
      		 
      		  bw.write(Data[i].toString().replace("[", "").replace("]", "").replace(" ", "") + "\r\n");      		  
      	  }
    	  bw.close();
          } catch (Exception e) {
          e.printStackTrace();
      }
	  long endTime = System.currentTimeMillis();
	  //System.out.println("Use " + (endTime - startTime)/1000.0 + "s.");
	  System.out.println("Done.");
	  System.out.println("Next step aproir.r");
  	} catch (IOException e) {
  		e.printStackTrace();
  	}
 }
}