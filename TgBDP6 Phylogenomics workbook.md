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

# Extended data set Phylogenic Tree

		 #Toxoplasma gondii ME49
          awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TGME49_311625.fa > Toxoplasma_gondii.fa
          awk '/^>/{print ">Toxoplasma_gondii_ME49" ++i; next}{print}' Toxoplasma_gondii.fa > header_Toxoplasma_gondii.fa
    
 	 # Hepatospora eriocheir strain canceri: A0H76_1719
  		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' A0H76_1719.fa > Hepatospora_eriocheir.fa
          awk '/^>/{print ">Hepatospora_eriocheir" ++i; next}{print}' Hepatospora_eriocheir.fa > header_Hepatospora_eriocheir.fa
         
	#Amphiamblys sp. WSBS2006: A8A55_1843
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' A8A55_1843.fa > Amphiamblys_sp.fa
          awk '/^>/{print ">Amphiamblys_sp" ++i; next}{print}' Amphiamblys_sp.fa > header_Amphiamblys_sp.fa
	
	#Nosema ceranae strain PA08_1199: AAJ76_200067555
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' AAJ76_200067555.fa > Nosema_ceranae.fa
          awk '/^>/{print ">Nosema_ceranae" ++i; next}{print}' Nosema_ceranae.fa > header_Nosema_ceranae.fa
	# Acanthamoeba castellanii str. Neff : ACA1_067440
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  ACA1_067440.fa > Acanthamoeba_castellanii.fa
          awk '/^>/{print ">Acanthamoeba_castellanii" ++i; next}{print}' Acanthamoeba_castellanii.fa > header_Acanthamoeba_castellanii.fa
	# Encephalitozoon cuniculi EC1: AEWD_051510
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}'  AEWD_051510.fa > Encephalitozoon_cuniculi.fa
          awk '/^>/{print ">Encephalitozoon_cuniculi" ++i; next}{print}' Encephalitozoon_cuniculi.fa > header_Encephalitozoon_cuniculi.fa
	# Babesia bigemina strain BOND: BBBOND_0100670
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BBBOND_0100670.fa > Babesia_bigemina.fa
          awk '/^>/{print ">Babesia_bigemina" ++i; next}{print}' Babesia_bigemina.fa > header_Babesia_bigemina.fa
	# Babesia bovis T2Bo: BBOV_IV004000
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BBOV_IV004000.fa > Babesia_bovis.fa
		  awk '/^>/{print ">Babesia_bovis" ++i; next}{print}' Babesia_bovis.fa > header_Babesia_bovis.fa
	# Theileria equi strain WA: BEWA_042410
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BEWA_042410.fa > Theileria_equi.fa
			  awk '/^>/{print ">Theileria_equi" ++i; next}{print}' Theileria_equi.fa > header_Theileria_equi.fa
	# Babesia microti strain RI: BMR1_01G00076
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BMR1_01G00076.fa > Babesia_microti.fa
			  awk '/^>/{print ">Babesia_microti" ++i; next}{print}' Babesia_microti.fa > header_Babesia_microti.fa
	# Babesia ovata strain Miyake: BOVATA_026070
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' BOVATA_026070.fa > Babesia_ovata.fa
				  awk '/^>/{print ">Babesia_ovata" ++i; next}{print}' Babesia_ovata.fa > header_Babesia_ovata.fa
	# Babesia divergens strain 1802A: Bdiv_029230
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Bdiv_029230.fa > Babesia_divergens.fa
					  awk '/^>/{print ">Babesia_divergens" ++i; next}{print}' Babesia_divergens.fa > header_Babesia_divergens.fa
	# Cryptosporidium ubiquitum isolate 39726 : cubi_02356 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cubi_02356.fa > Cryptosporidium_ubiquitum.fa
		 awk '/^>/{print ">Babesia_divergens" ++i; next}{print}' Babesia_divergens.fa > header_Babesia_divergens.fa
	# Cryptosporidium parvum Iowa II: cgd6_2320
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cgd6_2320.fa > Cryptosporidium_parvum.fa
			 awk '/^>/{print ">Cryptosporidium_parvum" ++i; next}{print}' Cryptosporidium_parvum.fa > header_Cryptosporidium_parvum.fa
	#Cryptosporidium andersoni isolate 30847: cand_006010
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' cand_006010.fa > Cryptosporidium_andersoni.fa
		awk '/^>/{print ">Cryptosporidium_andersoni" ++i; next}{print}' Cryptosporidium_andersoni.fa > header_Cryptosporidium_andersoni.fa
	# Vitrella brassicaformis CCMP3155: Vbra_11206
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Vbra_11206.fa > Vitrella_brassicaformis.fa
			awk '/^>/{print ">Vitrella_brassicaformis" ++i; next}{print}' Vitrella_brassicaformis.fa > header_Vitrella_brassicaformis.fa
	# Vittaforma corneae ATCC 50505: VICG_00309
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' VICG_00309.fa > Vittaforma_corneae.fa
		awk '/^>/{print ">Vittaforma_corneae" ++i; next}{print}' Vittaforma_corneae.fa > header_Vittaforma_corneae.fa
	# Theileria parva strain Muguga: TP01_1030
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TP01_1030.fa > Theileria_parva.fa
		awk '/^>/{print ">Theileria_parva" ++i; next}{print}' Theileria_parva.fa > header_Theileria_parva.fa
	# Theileria orientalis strain Shintoku: TOT_010000975
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TOT_010000975.fa > Theileria_orientalis.fa
		awk '/^>/{print ">Theileria_orientalis" ++i; next}{print}' Theileria_orientalis.fa > header_Theileria_orientalis.fa
	# Theileria annulata strain Ankara: TA16925
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' TA16925.fa > Theileria_annulata.fa
		awk '/^>/{print ">Theileria_annulata" ++i; next}{print}' Theileria_annulata.fa > header_Theileria_annulata.fa
	# Saprolegnia parasitica CBS 223.65: SPRG_12989
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SPRG_12989.fa > Saprolegnia_parasitica.fa
		awk '/^>/{print ">Saprolegnia_parasitica" ++i; next}{print}' Saprolegnia_parasitica.fa > header_Saprolegnia_parasitica.fa
	# Spraguea lophii 42_110: SLOPH_2073
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SLOPH_2073.fa > Spraguea_lophii.fa
		awk '/^>/{print ">Spraguea_lophii" ++i; next}{print}' Spraguea_lophii.fa > header_Spraguea_lophii.fa
	# Nematocida ausubeli ERTm2: NERG_01521
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NERG_01521.fa > Nematocida_ausubeli.fa
		awk '/^>/{print ">Nematocida_ausubeli" ++i; next}{print}' Nematocida_ausubeli.fa > header_Nematocida_ausubeli.fa
	# Cytauxzoon felis strain Winnie: CF001835
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CF001835.fa > Cytauxzoon_felis.fa
			awk '/^>/{print ">Cytauxzoon_felis" ++i; next}{print}' Cytauxzoon_felis.fa > header_Cytauxzoon_felis.fa
	#  Cryptosporidium hominis UdeA01: CHUDEA6_2320
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CHUDEA6_2320.fa > Cryptosporidium_hominis.fa
		awk '/^>/{print ">Cryptosporidium_hominis" ++i; next}{print}' Cryptosporidium_hominis.fa > header_Cryptosporidium_hominis.fa
	# Cryptosporidium muris RN66: CMU_007550
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CMU_007550.fa > Cryptosporidium_muris.fa
			awk '/^>/{print ">Cryptosporidium_muris" ++i; next}{print}' Cryptosporidium_muris.fa > header_Cryptosporidium_muris.fa
	# Cryptosporidium tyzzeri isolate UGA55: CTYZ_00001415
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CTYZ_00001415.fa > Cryptosporidium_tyzzeri.fa
			awk '/^>/{print ">Cryptosporidium_tyzzeri" ++i; next}{print}' Cryptosporidium_tyzzeri.fa > header_Cryptosporidium_tyzzeri.fa
	# Cryptosporidium meleagridis strain UKMEL1: CmeUKMEL1_12945
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CmeUKMEL1_12945.fa > Cryptosporidium_meleagridis.fa
	awk '/^>/{print ">Cryptosporidium_meleagridis" ++i; next}{print}' Cryptosporidium_meleagridis.fa > header_Cryptosporidium_meleagridis.fa
	# Chromera velia CCMP2878: Cvel_9760 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Cvel_9760.fa > Chromera_velia.fa
		awk '/^>/{print ">Chromera_velia" ++i; next}{print}' Chromera_velia.fa > header_Chromera_velia.fa
	# Mitosporidium daphniae UGP3: DI09_160p40
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' DI09_160p40.fa > Mitosporidium_daphniae.fa
		awk '/^>/{print ">Mitosporidium_daphniae" ++i; next}{print}' Mitosporidium_daphniae.fa > header_Mitosporidium_daphniae.fa
	# Enterocytozoon bieneusi H348: EBI_27476
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EBI_27476.fa > Enterocytozoon_bieneusi.fa
		awk '/^>/{print ">Enterocytozoon_bieneusi" ++i; next}{print}' Enterocytozoon_bieneusi.fa > header_Enterocytozoon_bieneusi.fa
	# Enterospora canceri strain GB1: ECANGB1_970
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ECANGB1_970.fa > Enterospora_canceri.fa
			awk '/^>/{print ">Enterospora_canceri" ++i; next}{print}' Enterospora_canceri.fa > header_Enterospora_canceri.fa
	# Encephalitozoon hellem ATCC 50504: EHEL_081170
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EHEL_081170.fa > Encephalitozoon_hellem.fa
			awk '/^>/{print ">Encephalitozoon_hellem" ++i; next}{print}' Encephalitozoon_hellem.fa > header_Encephalitozoon_hellem.fa
	# Entamoeba histolytica HM-1:IMSS: EHI_146100
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EHI_146100.fa > Entamoeba_histolytica.fa
				awk '/^>/{print ">Entamoeba_histolytica" ++i; next}{print}' Entamoeba_histolytica.fa > header_Entamoeba_histolytica.fa
	#  Enterocytozoon hepatopenaei strain TH1: EHP00_2617
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EHP00_2617.fa > Enterocytozoon_hepatopenaei.fa
				awk '/^>/{print ">Enterocytozoon_hepatopenaei" ++i; next}{print}' Enterocytozoon_hepatopenaei.fa > header_Enterocytozoon_hepatopenaei.fa
	# Entamoeba moshkovskii Laredo: EMO_100830 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EMO_100830.fa > Entamoeba_moshkovskii.fa
					awk '/^>/{print ">Entamoeba_moshkovskii" ++i; next}{print}' Entamoeba_moshkovskii.fa > header_Entamoeba_moshkovskii.fa
	# Mus musculus C57BL6J: ENSMUSG0000002291415
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENSMUSG0000002291415.fa > Mus_musculus.fa
					awk '/^>/{print ">Mus_musculus" ++i; next}{print}' Mus_musculus.fa > header_Mus_musculus.fa
	# Entamoeba nuttalli P19: ENU1_033640 
	awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' ENU1_033640.fa > Entamoeba_nuttalli.fa
					awk '/^>/{print ">Entamoeba_nuttalli" ++i; next}{print}' Entamoeba_nuttalli.fa > header_Entamoeba_nuttalli.fa
	# Nematocida parisii : NEPG_00799 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NEPG_00799.fa > Nematocida_parisii.fa
			awk '/^>/{print ">Nematocida_parisii" ++i; next}{print}' Nematocida_parisii.fa > header_Nematocida_parisii.fa
	# Nematocida displodere strain JUm2807: NEDG_00342
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NEDG_00342.fa > Nematocida_displodere.fa
				awk '/^>/{print ">Nematocida_displodere" ++i; next}{print}' Nematocida_displodere.fa > header_Nematocida_displodere.fa
	# Nosema bombycis CQ1: NBO_74g0002 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NBO_74g0002.fa > Nosema_bombycis.fa
		awk '/^>/{print ">Nosema_bombycis" ++i; next}{print}' Nosema_bombycis.fa > header_Nosema_bombycis.fa
	# Naegleria gruberi strain NEG-M: NAEGRDRAFT_78346 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NAEGRDRAFT_78346.fa > Naegleria_gruberi.fa
		awk '/^>/{print ">Naegleria_gruberi" ++i; next}{print}' Naegleria_gruberi.fa > header_Naegleria_gruberi.fa
	# Monocercomonoides exilis PA203: MONOS_15088 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' MONOS_15088.fa > Monocercomonoides_exilis.fa
			awk '/^>/{print ">Monocercomonoides_exilis" ++i; next}{print}' Monocercomonoides_exilis.fa > header_Monocercomonoides_exilis.fa
	# Encephalitozoon romaleae SJ-2008: EROM_081190 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' EROM_081190.fa > Encephalitozoon_romaleae.fa
			awk '/^>/{print ">Encephalitozoon_romaleae" ++i; next}{print}' Encephalitozoon_romaleae.fa > header_Encephalitozoon_romaleae.fa
	# Encephalitozoon intestinalis ATCC 50506: Eint_081180 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' Eint_081180.fa > Encephalitozoon_intestinalis.fa
		awk '/^>/{print ">Encephalitozoon_intestinalis" ++i; next}{print}' Encephalitozoon_intestinalis.fa > header_Encephalitozoon_intestinalis.fa
	# Aphanomyces astaci strain APO3: H257_13465 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' H257_13465.fa > Aphanomyces_astaci.fa
		awk '/^>/{print ">Aphanomyces_astaci" ++i; next}{print}' Aphanomyces_astaci.fa > header_Aphanomyces_astaci.fa
	# Aphanomyces invadans NJM9701: H310_06274 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' H310_06274.fa > Aphanomyces_invadans.fa
		awk '/^>/{print ">Aphanomyces_invadans" ++i; next}{print}' Aphanomyces_invadans.fa > header_Aphanomyces_invadans.fa
	# Anncaliia algerae PRA109: H311_04671 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' H311_04671.fa > Anncaliia_algerae.fa
		awk '/^>/{print ">Anncaliia_algerae" ++i; next}{print}' Anncaliia_algerae.fa > header_Anncaliia_algerae.fa
	# Ordospora colligata OC4: M896_091190 
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' M896_091190.fa > Ordospora_colligata.fa
		awk '/^>/{print ">Ordospora_colligata" ++i; next}{print}' Ordospora_colligata.fa > header_Ordospora_colligata.fa
	#Hammondia hammondi strain H.H.34: HHA_311625
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' HHA_311625.fa > Hammondia_hammondi.fa
		awk '/^>/{print "> Hammondia_hammondi " ++i; next}{print}' Hammondia_hammondi.fa > header_Hammondia_hammondi.fa
	# Cyclospora_cayetanensis: LOC34618198
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' LOC34618198.fa > Cyclospora_cayetanensis.fa
		awk '/^>/{print ">Cyclospora_cayetanensis " ++i; next}{print}' Cyclospora_cayetanensis.fa > header_Cyclospora_cayetanensis.fa
	# Sarcocystis_neurona: SN3_00801045  
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' SN3_00801045.fa > Sarcocystis_neurona.fa
		awk '/^>/{print ">Sarcocystis_neurona " ++i; next}{print}' Sarcocystis_neurona.fa > header_Sarcocystis_neurona.fa
	# Cystoisospora_suis: CSUI_000595  
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' CSUI_000595.fa > Cystoisospora_suis.fa
		awk '/^>/{print ">Cystoisospora_suis" ++i; next}{print}' Cystoisospora_suis.fa > header_Cystoisospora_suis.fa
	# Neospora_caninum: NCLIV_055410  
		awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' NCLIV_055410.fa > Neospora_caninum.fa
		awk '/^>/{print ">Neospora_caninum " ++i; next}{print}' Neospora_caninum.fa > header_Neospora_caninum.fa

		cat header_Toxoplasma_gondii.fa header_Hepatospora_eriocheir.fa header_Amphiamblys_sp.fa header_Nosema_ceranae.fa header_Acanthamoeba_castellanii.fa header_Encephalitozoon_cuniculi.fa header_Babesia_bigemina.fa header_Babesia_bovis.fa header_Theileria_equi.fa header_Babesia_microti.fa header_Babesia_ovata.fa header_Babesia_divergens.fa header_Cryptosporidium_ubiquitum.fa header_Cryptosporidium_parvum.fa header_Cryptosporidium_andersoni.fa header_Vitrella_brassicaformis.fa header_Vittaforma_corneae.fa header_Theileria_parva.fa header_Theileria_orientalis.fa header_Theileria_annulata.fa header_Saprolegnia_parasitica.fa header_Spraguea_lophii.fa header_Nematocida_ausubeli.fa header_Cytauxzoon_felis.fa header_Cryptosporidium_hominis.fa header_Cryptosporidium_muris.fa header_Cryptosporidium_tyzzeri.fa header_Cryptosporidium_meleagridis.fa header_Chromera_velia.fa header_Mitosporidium_daphniae.fa header_Enterocytozoon_bieneusi.fa header_Enterospora_canceri.fa header_Encephalitozoon_hellem.fa header_Entamoeba_histolytica.fa header_Enterocytozoon_hepatopenaei.fa header_Entamoeba_moshkovskii.fa header_Mus_musculus.fa header_Entamoeba_nuttalli.fa header_Nematocida_parisii.fa header_Nematocida_displodere.fa header_Nosema_bombycis.fa header_Naegleria_gruberi.fa header_Monocercomonoides_exilis.fa header_Encephalitozoon_romaleae.fa header_Encephalitozoon_intestinalis.fa header_Aphanomyces_astaci.fa header_Aphanomyces_invadans.fa header_Anncaliia_algerae.fa header_Ordospora_colligata.fa header_Hammondia_hammondi.fa header_Cyclospora_cayetanensis.fa header_Sarcocystis_neurona.fa header_Cystoisospora_suis.fa header_Neospora_caninum.fa > TGBDP6_Phylogenetics_Trees_4.fa

		mafft --auto TGBDP6_Phylogenetics_Trees_4.fa > Output_TGBDP6_Phylogenetics_Tree_4.fa

	iqtree -s Output_TGBDP6_Phylogenetics_Tree_4.fa -m LG -bb 1000 -pre output
	
	
	
	(Toxoplasma_gondii_ME491:0.057030,(((((((((((Hepatospora_eriocheir1:1.301217,((Nosema_ceranae1:0.680232,(((Encephalitozoon_hellem1:0.094674,Encephalitozoon_romaleae1:0.095216)100:0.102742,Encephalitozoon_intestinalis1:0.142409)100:0.203449,Ordospora_colligata1:0.453473)100:0.269436)100:0.413591,Spraguea_lophii1:0.796900)96:0.101085)98:0.183650,(Nematocida_ausubeli1:0.606900,Nematocida_displodere1:0.483357)100:0.633691)100:0.569604,(Encephalitozoon_cuniculi1:1.398377,Nematocida_parisii1:1.205662)65:0.112909)55:0.204505,((((Vittaforma_corneae1:0.476167,((Enterocytozoon_bieneusi1:0.517550,Enterocytozoon_hepatopenaei1:0.578263)100:0.164537,Enterospora_canceri1:0.532862)96:0.204160)86:0.155092,Nosema_bombycis1:0.416065)100:0.929232,(Naegleria_gruberi1:0.781587,Anncaliia_algerae1:1.004465)49:0.111548)92:0.282221,(((Entamoeba_histolytica1:1.223560,Entamoeba_nuttalli1:0.814689)97:0.161915,Entamoeba_moshkovskii1:1.734440)92:0.146152,Monocercomonoides_exilis1:1.152767)71:0.177524)34:0.087964)32:0.096195,(((Amphiamblys_sp1:1.238953,(Acanthamoeba_castellanii1:0.725624,Mus_musculus1:0.887513)97:0.109057)90:0.131030,(Saprolegnia_parasitica1:0.405943,(Aphanomyces_astaci1:0.102801,Aphanomyces_invadans1:0.140940)100:0.351813)100:0.502824)99:0.293034,Mitosporidium_daphniae1:1.350634)86:0.194441)48:0.166573,(Vitrella_brassicaformis1:0.924986,Chromera_velia1:1.068773)100:0.200146)46:0.119848,((((((Babesia_bigemina1:0.122408,Babesia_ovata1:0.108050)100:0.324577,Babesia_divergens1:0.472340)100:0.135779,Babesia_bovis1:0.580375)100:0.317534,((Theileria_equi1:0.378591,Cytauxzoon_felis1:0.529194)100:0.231662,((Theileria_parva1:0.151126,Theileria_annulata1:0.099520)100:0.580741,Theileria_orientalis1:0.500975)100:0.370983)100:0.220981)100:0.302498,Babesia_microti1:1.116078)100:0.363992,((Cryptosporidium_ubiquitum:0.103873,(((Cryptosporidium_parvum1:0.009774,Cryptosporidium_tyzzeri1:0.010627)100:0.004118,Cryptosporidium_hominis1:0.010408)100:0.020249,Cryptosporidium_meleagridis1:0.042730)100:0.112764)100:0.412869,(Cryptosporidium_andersoni1:0.017024,Cryptosporidium_muris1:0.014860)100:0.562090)100:0.761307)79:0.176888)100:0.318718,Cyclospora_cayetanensis:0.822854)100:0.156207,Sarcocystis_neurona:0.691418)100:0.101290,Cystoisospora_suis1:0.725144)100:0.470932,Neospora_caninum:0.225331)100:0.180525,Hammondia_hammondi:0.048571);
	cat output.contree
