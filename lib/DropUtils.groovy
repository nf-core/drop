import groovy.yaml.YamlSlurper
import java.nio.file.Path

public class Utils {
    public static Map readYamlFile(Path yaml) {
        return new YamlSlurper().parse(yaml.text)
    }
}
