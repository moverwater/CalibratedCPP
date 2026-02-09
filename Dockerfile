# Use an official Maven image with Java 21
FROM maven:3.9.6-eclipse-temurin-21

WORKDIR /app

# Copy the project files
COPY calibratedcpp-beast/ .

# 1. REMOVE the broken parent block
RUN sed -i '/<parent>/,/<\/parent>/d' pom.xml

# 2. INJECT the missing groupId and Version
RUN sed -i '/<artifactId>calibratedcpp-beast<\/artifactId>/i <groupId>io.github.linguaphylo<\/groupId>' pom.xml && \
    sed -i 's/<version>${revision}<\/version>/<version>0.0.1<\/version>/' pom.xml

# 3. INJECT JUnit 5 Dependency (The Fix)
# We find the <dependencies> tag and append the JUnit XML block right after it.
RUN sed -i '/<dependencies>/a \
    <dependency> \
        <groupId>org.junit.jupiter</groupId> \
        <artifactId>junit-jupiter</artifactId> \
        <version>5.10.2</version> \
        <scope>test</scope> \
    </dependency>' pom.xml

# 4. Install local JARs
RUN mvn install:install-file -Dfile=lib/BEAST.base-2.7.8.jar -DgroupId=beast2 -DartifactId=beast-base -Dversion=2.7.8 -Dpackaging=jar && \
    mvn install:install-file -Dfile=lib/BEAST.app-2.7.8.jar -DgroupId=beast2 -DartifactId=beast-app -Dversion=2.7.8 -Dpackaging=jar && \
    mvn install:install-file -Dfile=lib/launcher-2.7.8.jar -DgroupId=beast2 -DartifactId=beast-launcher -Dversion=2.7.8 -Dpackaging=jar

# 5. Run the tests
RUN mvn clean test