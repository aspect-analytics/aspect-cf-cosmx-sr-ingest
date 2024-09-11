#!groovy
// Pipeline configration for CI/CD
// For this projecte the build results can be found at https://jenkins.aspect.aspect-analytics.ninja/job/github/job/python_project_template/

@Library(['jpipe', 'aspect'])

String appName = 'aspect-py-project-template'

// Pipeline configration
// For more information see: https://docs.aspect.aspect-analytics.ninja/pipeline/jenkinsfile/
Map config = [
  services: [
    [
      // Source
      imageName    : appName,  // Name of the container imag
      dockerfile   : './infra/Dockerfile',
      // Publish Python Package
      // Packages published at https://nexus.aspect.aspect-analytics.ninja/#browse/browse:pypi-aspect-analytics
      pythonPublish: true,  // Publish package to internal repository, set to true if you want to publish to our internal repository
      pythonDistPath: '/app/repository/dist',  // Expect this to be provided by the container build
      // Publish Docker Image
      imagePublish: false,  // No need to publish images for packages
    ]
  ],
]

aspectProject config
