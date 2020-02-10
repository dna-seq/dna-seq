version development

workflow hello {
    input {
        String name
    }

    call greeting {
        input:
            name = name
    }

    output {
        String out = greeting.out
    }

}

task greeting {
    input {
        String name
    }

    command {
        echo "Hello ~{name}!"
    }

    output {
        String out = read_string(stdout())
    }
}
