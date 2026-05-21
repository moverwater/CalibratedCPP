# Force standard Intel architecture so JavaFX 17 can be downloaded via Rosetta 2
FROM maven:3.9.16-eclipse-temurin-25

# Set working directory
WORKDIR /app

# Copy ONLY the beast module, ignoring the root project that requires LPhyBeast
COPY calibratedcpp-beast/ .

# Install a proper XML manipulation tool
RUN apt-get update && apt-get install -y xmlstarlet

# Safely edit the XML tree: Delete the parent block and add the groupId
RUN xmlstarlet ed --inplace \
    -d "/project/parent" \
    -s "/project" -t elem -n "groupId" -v "io.github.linguaphylo" \
    pom.xml

# Run the tests specifically for the beast module
RUN mvn clean test -Drevision=0.0.1