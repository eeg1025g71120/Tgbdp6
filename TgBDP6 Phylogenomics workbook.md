# Analyzing TgBDP6 through BLAST and Phylogenomic Analyses to Access Evolutionary Relationships Amongst Apicomplexians
This is a workbook detailing the methods and steps utilized in this bioinformatic project to access if the TgBDP6 sequence is conserved across the Apicomplexan phylum.
# Background

# Materials
   * PC: Lenova Yoga 730
   * Coding Terminal: PuTTy
   * Preformance Cluster : Ron
   * Programs Utalized:
       * MAFFT
       * iqTree
       * FigTree v1.4.4 ( https://github.com/rambaut/figtree/releases)
    * Databases Utalized:
      *  VEuPathDB (https://veupathdb.org/veupathdb/app)
      *   UniProtKB (https://www.uniprot.org/)
# Results from BLAST fr EUPATHdb
* 	Hammondia hammondi strain H.H.34: HHA_311625
*   Besnoitia besnoiti strain Bb-Ger1: BESB_072950
*	Cystoisospora suis strain Wien I: CSUI_000595
* Sarcocystis neurona SN3: SN3_00801045
* Cyclospora cayetanensis isolate NF1_C8: LOC34618198
* Eimeria falciformis Bayer Haberkorn 1970: EfaB_MINUS_25458
* Eimeria maxima Weybridge: EMWEY_00047460
* Eimeria acervulina Houghton: EAH_00038960
* Cryptosporidium muris RN66: CMU_007550
* Cryptosporidium andersoni isolate 30847: cand_006010
* Babesia bigemina strain BOND: BBBOND_0100670
* Theileria orientalis strain Shintoku: TOT_010000975
* Babesia bovis T2Bo: BBOV_IV004000
* Theileria equi strain WA BEWA_042410
* Cryptosporidium hominis UdeA01 CHUDEA6_2320
* Cryptosporidium parvum IOWA-ATCC  CPATCC_0013840
* Cryptosporidium tyzzeri isolate UGA55: CTYZ_00001415
* Cryptosporidium ubiquitum isolate 39726:cubi_02356
* Cryptosporidium meleagridis strain UKMEL1: CmeUKMEL1_12945
* Babesia divergens strain 1802A: Bdiv_029230
* Eimeria necatrix Houghton ENH_00053350
* Cytauxzoon felis strain Winnie CF001835
# Setting up for the Phylogenomic Analysis:
First I created a new directory to begin my analysis in:
	
	mkdir TgBDP6
		
Next, I moved into the newly created directed to begin my analysis:
	
	cd TgBDP6

# Phylogenomics Analysis of TGBDP4 predicted protein sequence within Aplicomplexians 

      # Phylogenetic
          #Toxoplasma gondii ME49
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TGME49_311625.fa > Toxoplasma_gondii.fa
              awk '/^>/{print ">Toxoplasma_gondii_ME49" ++i; next}{print}' Toxoplasma_gondii.fa > header_Toxoplasma_gondii.fa
          # Hammondia hammondi strain H.H.34: HHA_311625
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' HHA_311625.fa > Hammondia_hammondi.fa
              awk '/^>/{print "> Hammondia_hammondi_HH34" ++i; next}{print}' Hammondia_hammondi.fa > header_Hammondia_hammondi.fa
          # Besnoitia besnoiti strain Bb-Ger1: BESB_072950
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BESB_072950.fa > Besnoitia_besnoiti.fa
              awk '/^>/{print "> Besnoitia_besnoiti_Bb-Ger1" ++i; next}{print}' Besnoitia_besnoiti.fa > header_Besnoitia_besnoiti.fa
           # Cystoisospora suis strain Wien I: CSUI_000595
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CSUI_000595 .fa > Cystoisospora_suis.fa
              awk '/^>/{print "> Cystoisospora_suis" ++i; next}{print}' Cystoisospora_suis.fa > header_Cystoisospora_suis.fa
          # Sarcocystis neurona SN3: SN3_00801045
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SN3_00801045.fa > Sarcocystis_neurona.fa
              awk '/^>/{print "> Sarcocystis_neurona" ++i; next}{print}' Sarcocystis_neurona.fa > header_Sarcocystis_neurona.fa
          # Cyclospora cayetanensis isolate NF1_C8: LOC34618198
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' LOC34618198.fa > Cyclospora_cayetanensis.fa
              awk '/^>/{print "> Cyclospora_cayetanensis" ++i; next}{print}' Cyclospora_cayetanensis.fa > header_Cyclospora_cayetanensis.fa
          # Eimeria falciformis Bayer Haberkorn 1970: EfaB_MINUS_25458
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EfaB_MINUS_25458.fa> Eimeria_falciformis.fa
              awk '/^>/{print "> Eimeria_falciformis" ++i; next}{print}' Eimeria_falciformis.fa > header_Eimeria_falciformis.fa
          # Eimeria maxima Weybridge: EMWEY_00047460
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EMWEY_00047460.fa > Eimeria_maxima.fa
              awk '/^>/{print "> Eimeria_maxima" ++i; next}{print}' Eimeria_maxima.fa > header_Eimeria_maxima.fa
          # Eimeria acervulina Houghton: EAH_00038960
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EAH_00038960.fa > Eimeria_acervulina.fa
              awk '/^>/{print "> Eimeria_acervulina" ++i; next}{print}' Eimeria_acervulina.fa > header_Eimeria_acervulina.fa
          # Cryptosporidium muris RN66: CMU_007550
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CMU_007550.fa > Cryptosporidium_muris.fa
              awk '/^>/{print "> Cryptosporidium_muris" ++i; next}{print}' Cryptosporidium_muris.fa > header_Cryptosporidium_muris.fa
          # Cryptosporidium andersoni isolate 30847: cand_006010
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cand_006010.fa > Cryptosporidium_andersoni.fa
		awk '/^>/{print "> Cryptosporidium_andersoni" ++i; next}{print}' Cryptosporidium_andersoni.fa > header_Cryptosporidium_andersoni.fa
          # Babesia bigemina strain BOND: BBBOND_0100670
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'BBBOND_0100670.fa > Babesia_bigemina.fa
		awk '/^>/{print "> Babesia_bigemina" ++i; next}{print}' Babesia_bigemina.fa > header_Babesia_bigemina.fa
          # Theileria orientalis strain Shintoku: TOT_010000975
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TOT_010000975.fa > Theileria_orientalis.fa
		awk '/^>/{print "> Theileria_orientalis" ++i; next}{print}' Theileria_orientalis.fa > header_Theileria_orientalis.fa
          # Babesia bovis T2Bo: BBOV_IV004000
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BBOV_IV004000.fa > Babesia_bovis.fa
		awk '/^>/{print "> Babesia_bovis" ++i; next}{print}' Babesia_bovis.fa > header_Babesia_bovis.fa
          # Theileria equi strain WA BEWA_042410
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BEWA_042410.fa > Theileria_equi.fa
		awk '/^>/{print "> Theileria_equi" ++i; next}{print}' Theileria_equi.fa > header_Theileria_equi.fa
          # Cryptosporidium hominis UdeA01 CHUDEA6_2320
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CHUDEA6_2320.fa > Cryptosporidium_hominis.fa
		awk '/^>/{print "> Cryptosporidium_hominis" ++i; next}{print}' Cryptosporidium_hominis.fa > header_Cryptosporidium_hominis.fa
          # Cryptosporidium parvum IOWA-ATCC  CPATCC_0013840
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CPATCC_0013840.fa > Cryptosporidium_parvum.fa
		awk '/^>/{print "> Cryptosporidium_parvum" ++i; next}{print}' Cryptosporidium_parvum.fa > header_Cryptosporidium_parvum.fa
          # Cryptosporidium tyzzeri isolate UGA55: CTYZ_00001415
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CTYZ_00001415.fa > Cryptosporidium_tyzzeri.fa
		awk '/^>/{print "> Cryptosporidium_tyzzeri " ++i; next}{print}' Cryptosporidium_tyzzeri.fa > header_Cryptosporidium_tyzzeri.fa
          # Cryptosporidium ubiquitum isolate 39726:cubi_02356
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cubi_02356.fa > Cryptosporidium_ubiquitum.fa
		awk '/^>/{print "> Cryptosporidium_ubiquitum " ++i; next}{print}' Cryptosporidium_ubiquitum.fa > header_Cryptosporidium_ubiquitum.fa
          # Cryptosporidium meleagridis strain UKMEL1: CmeUKMEL1_12945
	  	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CmeUKMEL1_12945.fa > Cryptosporidium_meleagridis.fa
		awk '/^>/{print "> Cryptosporidium_meleagridis" ++i; next}{print}' Cryptosporidium_meleagridis.fa > header_Cryptosporidium_meleagridis.fa
          # Babesia divergens strain 1802A: Bdiv_029230
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Bdiv_029230.fa > Babesia_divergens.fa
              awk '/^>/{print "> Babesia_divergens" ++i; next}{print}' Babesia_divergens.fa > header_Babesia_divergens.fa
          # Eimeria necatrix Houghton ENH_00053350
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENH_00053350.fa > Eimeria_necatrix.fa
              awk '/^>/{print "> Eimeria_necatrix" ++i; next}{print}' Eimeria_necatrix.fa > header_Eimeria_necatrix.fa
          # Cytauxzoon felis strain Winnie CF001835
              awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CF001835.fa > Cytauxzoon_felis.fa
              awk '/^>/{print "> Cytauxzoon_felis " ++i; next}{print}' Cytauxzoon_felis.fa > header_Cytauxzoon_felis.fa

        cat header_Toxoplasma_gondii.fa header_Hammondia_hammondi.fa header_Besnoitia_besnoiti.fa header_Cystoisospora_suis.fa header_Sarcocystis_neurona.fa header_Cyclospora_cayetanensis.fa header_Eimeria_falciformis.fa header_Eimeria_maxima.fa header_Eimeria_acervulina.fa header_Cryptosporidium_muris.fa header_Cryptosporidium_andersoni.fa header_Babesia_bigemina.fa header_Theileria_orientalis.fa header_Babesia_bovis.fa header_Theileria_equi.fa header_Cryptosporidium_hominis.fa header_Cryptosporidium_parvum.fa header_Cryptosporidium_tyzzeri.fa header_Cryptosporidium_ubiquitum.fa header_Cryptosporidium_meleagridis.fa header_Babesia_divergens.fa header_Eimeria_necatrix.fa header_Cytauxzoon_felis.fa > TGBDP6_Phylogenetics_Tree.fa

	mafft --auto TGBDP6_Phylogenetics_Tree.fa > Output_TGBDP6_Phylogenetics_Tree.fa
	iqtree -s Output_TGBDP6_Phylogenetics_Tree.fa -m LG -bb 1000 -pre output
	cat output.contree
	
	
	
	
	(Toxoplasma_gondii_ME491:0.046906,Hammondia_hammondi_HH341:0.058573,(Besnoitia_besnoiti_Bb-Ger11:0.347050,(Sarcocystis_neurona1:0.566032,(((Cyclospora_cayetanensis1:0.336473,Eimeria_falciformis1:0.247412)85:0.071468,((Eimeria_maxima1:0.310487,Eimeria_acervulina1:0.149560)87:0.038062,Eimeria_necatrix1:0.193124)90:0.059766)100:0.483250,(((Cryptosporidium_muris1:0.015800,Cryptosporidium_andersoni1:0.016037)100:0.517879,(((Cryptosporidium_hominis1:0.010400,(Cryptosporidium_parvum1:0.009772,Cryptosporidium_tyzzeri:0.010623)99:0.004128)100:0.020289,Cryptosporidium_meleagridis1:0.042716)100:0.108773,Cryptosporidium_ubiquitum:0.108001)100:0.449123)100:0.778738,(((Babesia_bigemina1:0.379394,Babesia_divergens1:0.461547)100:0.153661,Babesia_bovis1:0.555256)100:0.290630,(Theileria_orientalis1:0.832565,(Theileria_equi1:0.394856,Cytauxzoon_felis:0.518317)100:0.206829)100:0.254584)100:0.583447)100:0.537649)99:0.190426)100:0.430465)100:0.354310);

# Results from NCBI
