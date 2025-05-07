import groovy.yaml.YamlSlurper

public static class Utils {
    public static Map readYamlFile(Path yaml) {
        return new YamlSlurper().parse(yaml.text)
    }
}
