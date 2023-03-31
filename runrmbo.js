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
        {   "test"       : "librembo/bin/test_rembo.x",
            "rembo"      : "librembo/bin/test_rembo.x"
        },

    "grp_aliases" : {},

    "par_paths" : {},

    "files" : [], 

    "dir-special" : {},

    "links" : 
        ["input","ice_data","maps"],

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
