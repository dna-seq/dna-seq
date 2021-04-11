version development

workflow opencravat {
 input {
     File modules
 }
}

    task annotate {
        input {
            File modules
        }

        command {
            oc
        }

        runtime {
            docker: "mcr.microsoft.com/opencravat"
        }
    }