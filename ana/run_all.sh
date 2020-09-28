#/usr/bin/env bash

# Directory where the simulated data samples are stored
data_dir=./rootfile.ptgun.piK

# remove file
rm info_input.txt

# flag for each momentum 
mom_0_5GeV=1
mom_0_75GeV=1
mom_1GeV=1
mom_3GeV=1
mom_5GeV=1
mom_10GeV=1
mom_50GeV=1
mom_100GeV=1

# flag for pion / kaon
ana_pion_flag=1
ana_kaon_flag=1

####################### Pion ########################

if [ ${ana_pion_flag} = 1 ];
then
    # remove file
    rm table_pion.txt
    
    # Pion 0.5GeV
    if [ ${mom_0_5GeV} = 1 ];
    then
        echo ${data_dir}/pion_0.5GeV_85deg_2000event.root   pion_0.50GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Pion 0.75GeV
    if [ ${mom_0_75GeV} = 1 ];
    then
        echo ${data_dir}/pion_0.75GeV_85deg_2000event.root   pion_0.75GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Pion 1GeV
    if [ ${mom_1GeV} = 1 ];
    then
        echo ${data_dir}/pion_1GeV_85deg_2000event.root   pion_1.00GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Pion 3GeV
    if [ ${mom_3GeV} = 1 ];
    then    
        echo ${data_dir}/pion_3GeV_85deg_2000event.root   pion_3.00GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Pion 5GeV
    if [ ${mom_5GeV} = 1 ];
    then
        echo ${data_dir}/pion_5GeV_85deg_2000event.root   pion_5.00GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi;

    # Pion 10GeV
    if [ ${mom_10GeV} = 1 ];
    then    
        echo ${data_dir}/pion_10GeV_85deg_2000event.root   pion_10.0GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Pion 50GeV
    if [ ${mom_50GeV} = 1 ];
    then    
        echo ${data_dir}/pion_50GeV_85deg_2000event.root   pion_50.0GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Pion 100GeV
    if [ ${mom_100GeV} = 1 ];
    then    
        echo ${data_dir}/pion_100GeV_85deg_2000event.root   pion_100.0GeV  table_pion.txt > info_input.txt
        root -l -q calc_dedx.c
    fi
fi


####################### Kaon ##########################

if [ ${ana_kaon_flag} = 1 ];
then
    # remove file
    rm table_kaon.txt
    
    # Kaon 0.5GeV
    if [ ${mom_0_5GeV} = 1 ];
    then    
        echo ${data_dir}/kaon_0.5GeV_85deg_2000event.root   kaon_0.50GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Kaon 0.75GeV
    if [ ${mom_0_75GeV} = 1 ];
    then    
        echo ${data_dir}/kaon_0.75GeV_85deg_2000event.root   kaon_0.75GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Kaon 1GeV
    if [ ${mom_1GeV} = 1 ];
    then    
        echo ${data_dir}/kaon_1GeV_85deg_2000event.root   kaon_1.00GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Kaon 3GeV
    if [ ${mom_3GeV} = 1 ];
    then    
        echo ${data_dir}/kaon_3GeV_85deg_2000event.root   kaon_3.00GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Kaon 5GeV
    if [ ${mom_5GeV} = 1 ];
    then    
        echo ${data_dir}/kaon_5GeV_85deg_2000event.root   kaon_5.00GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Kaon 10GeV
    if [ ${mom_10GeV} = 1 ];
    then    
        echo ${data_dir}/kaon_10GeV_85deg_2000event.root   kaon_10.0GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Kaon 50GeV
    if [ ${mom_50GeV} = 1 ];
    then        
        echo ${data_dir}/kaon_50GeV_85deg_2000event.root   kaon_50.0GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi

    # Kaon 100GeV
    if [ ${mom_100GeV} = 1 ];
    then        
        echo ${data_dir}/kaon_100GeV_85deg_2000event.root   kaon_100.0GeV  table_kaon.txt > info_input.txt
        root -l -q calc_dedx.c
    fi
    
fi



### making dE/dx distribution plots etc.
make_plot_flag=1

if [ ${make_plot_flag} = 1 ];
then
    root -l -q plot_dedx.c
fi


exit 0
