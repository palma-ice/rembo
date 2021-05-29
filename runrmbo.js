{
    "defaults" :
    {
        "jobname"       : "rembo",
        "email"         : "USER@pik-potsdam.de",
        "group"         : "anthroia", 
        "omp"           : 0,
        "wall"          : 24, 
        "qos"           : "priority",
        "partition"     : "haswell",
        "job_template"  : "config/pik_submit_slurm"
    },

    "exe_aliases" :
        {   "test"       : "librembo/bin/rembo_test.x",
            "rembo"      : "librembo/bin/rembo_test.x"
        },

    "grp_aliases" : {},

    "par_paths" : {},

    "files" : [], 

    "dir-special" : {},

    "links" : 
        ["input","ice_data","maps"],

    "const_paths" :
        {   "None"  : "None"
        },

    "const_path_default" : "par/rembo_const_Earth.nml",

    "job_queues" :
        {   "priority" :
            {   "wall" : 24  },
            "short" :
            {   "wall" : 24  },
            "medium" :
            {   "wall" : 168 },
            "long" :
            {   "wall" : 720 }
        }
}
